=head1 NAME

LGTSeek - Find Lateral Gene Transfer in sequencing data

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run computes and process the output of data for purposes
of finding putative lateral gene transfer.

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package LGTSeek;
use strict;
use version;
use File::Basename;
# Dependencies
use GiTaxon;
use LGTBestBlast;
use LGTFinder;
use LGTbwa;
use LGTsam2lca;
$| = 1;

=head2 new

 Title   : new
 Usage   : my $lgtseek = LGTSeek->new({fastq => $fastq,...})
 Function: Creates a new LGTSeek object.
 Returns : An instance of LGTSeek
 Args    : A hash containing potentially several config options:

           fastq - Path fastq files.
           host_fasta - One or more fasta files to use a host references
           donor_fasta - One or more fasta files to use as donor references
           prinseq_bin - Path to prinseq perl script
           bwa_bin - Path to bwa binary
           aspera_path - Path to Aspera ascp binary
           taxon_dir - Directory where taxon information lives
           taxon_idx_dir - Directory where taxon indices can be created (or where they are)
           taxon_host - Hostname of taxon mongo database

=cut

sub new {
  my($class,$args) = @_;

  my $self = $args;

  bless $self;
  return $self;
}

=head2 getGiTaxon

 Title   : getGiTaxon
 Usage   : my $gi2tax = $lgtseek->getGiTaxon({'host' => 'foobar.com'...});
 Function: Retrieve a GiTaxon object ready to assign taxonomic information
 Returns : A GiTaxon object
 Args    : A hash options to pass to GiTaxon. These could have been
           passed in globally. They can also be overridden here.

           taxon_dir - Directory where taxon information lives
           taxon_idx_dir - Directory where taxon indices can be created (or where they are)
           taxon_host - Hostname of taxon mongo database

=cut

sub getGiTaxon {
    my($self,$config) = @_;

    # If we already have a gitaxon object we'll just return it.
    if(!$self->{gitaxon}) {
        # Apply any config options that came over.
        if($config) {
            map {$self->{$_} = $config->{$_}} keys %$config;
        }

        # Create the object.
        $self->{gitaxon} = GiTaxon->new(
            {'taxon_dir' => $self->{taxon_dir},
             'chunk_size' => 10000,
             'idx_dir' => $self->{taxon_idx_dir},
             'host' => $self->{taxon_host},
             'type' => 'nucleotide'});        
    }
    
    return $self->{gitaxon};

}

=head2 prinseqFilterBam

 Title   : prinseqFilterBam
 Usage   : my $filteredBam = $LGTSeek->prinseqFilterBam({'bam_file' => '/path/to/file.bam'...})
 Function: Prinseq filter a bam file
 Returns : Path to the filtered bam file
 Args    : A bam_file and optionally the path to the prinseq_bin if this value
           was not already passed in.

=cut

sub prinseqFilterBam {
    my($self,$config) = @_;
    
    # Override if it is provided
    $self->{prinseq_bin} = $config->{prinseq_bin} ? $config->{prinseq_bin} :$self->{prinseq_bin};
    $self->{samtools_bin} = $self->{samtools_bin} ? $self->{samtools_bin} : 'samtools';
    if($config->{output_dir}) {
        $self->_run_cmd("mkdir -p $config->{output_dir}");
    }
    
    if(!$self->{ergatis_bin} || !$self->{prinseq_bin}) {
        die "Must provide an ergatis_bin and prinseq_bin parameter to run prinseq filtering $self->{prinseq_bin} $self->{ergatis_bin}\n";
    }

    my $retval;
    if($self->{paired_end}) {
        $retval = $self->_prinseqFilterPaired($config->{bam_file},$config->{output_dir});
    }
    else {
        die "Single end is currently not implemented\n";
    }

    return $retval;
}

=head2 _prinseqFilterPaired

 Title   : _prinseqFilterPaired
 Usage   : *PRIVATE*
 Function: Prinseq filter a bam paired end file
 Returns : Path to the filtered bam file
 Args    : A bam file to filter

=cut

sub _prinseqFilterPaired {
    my($self,$bam_file,$output_dir) = @_;

    my($name,$path,$suff) = fileparse($bam_file,'.bam');

    $output_dir = $output_dir ? $output_dir : $path;

    my $bin = $self->{bin_dir};
    my $prinseq_bin = $self->{prinseq_bin};

    my $samtools = $self->{samtools_bin};

    # Generate concatenated fastq files for prinseq derep filtering
    my $cmd = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$bam_file --fastq=1 --combine_mates=1 --output_file=$output_dir/$name\_combined.fastq";
    print STDERR "$cmd\n";
    $self->_run_cmd($cmd);
    # Run prinseq for dereplication
    my $cmd = "perl $prinseq_bin --fastq=$output_dir/$name\_combined.fastq --out_good=$output_dir/$name\_derep_good --out_bad=$output_dir/$name\_derep_bad -derep 14";
    print STDERR "$cmd\n";
    $self->_run_cmd($cmd);

    # Pull out bad ids
    my $cmd = "perl -e 'while(<>){s/\@//;print;<>;<>;<>;}' $output_dir/$name\_derep_bad.fastq > $output_dir/$name\_derep_bad_ids.out";
    print STDERR "$cmd\n";
    $self->_run_cmd($cmd);

    # Generate single-read fastq for low complexity filtering
    my $cmd = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$bam_file --fastq=1 --combine_mates=0 --paired=1 --output_file=$output_dir/$name.fastq";
    $self->_run_cmd($cmd);

    # Run prinseq for low complexity filtering
    my $cmd = "perl $prinseq_bin --fastq=$output_dir/$name\_1.fastq --out_good=$output_dir/$name\_lc_1_good --out_bad=$output_dir/$name\_lc_1_bad -lc_method dust -lc_threshold 7";
    $self->_run_cmd($cmd);

    if( -e "$output_dir/$name\_lc_1_bad.fastq") {
        # Pull out bad ids
        my $cmd = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $output_dir/$name\_lc_1_bad.fastq > $output_dir/$name\_lc_1_bad_ids.out";
        $self->_run_cmd($cmd);
    }
    else {
        print STDERR "Didn't find any low complexity sequences in read 1\n";
        $self->_run_cmd("touch $output_dir/$name\_lc_1_bad_ids.out");
    }

    # Run prinseq for low complexity filtering
    my $cmd = "perl $prinseq_bin --fastq=$output_dir/$name\_2.fastq --out_good=$output_dir/$name\_lc_2_good --out_bad=$output_dir/$name\_lc_2_bad -lc_method dust -lc_threshold 7";
    $self->_run_cmd($cmd);

    # Pull out bad ids
    if( -e "$output_dir/$name\_lc_2_bad.fastq") {
        my $cmd = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $output_dir/$name\_lc_2_bad.fastq > $output_dir/$name\_lc_2_bad_ids.out";
        $self->_run_cmd($cmd);
    }
    else {
        print STDERR "Didn't find any low complexity sequences in read 2\n";
        $self->_run_cmd("touch $output_dir/$name\_lc_2_bad_ids.out");
    }

    # Merge bad ids from derep and lc filtering
    my $cmd = "cat $output_dir/$name\_derep_bad_ids.out $output_dir/$name\_lc_1_bad_ids.out $output_dir/$name\_lc_2_bad_ids.out | sort -u > $output_dir/$name\_bad_ids.out";
   $self->_run_cmd($cmd);

    # Filter sam file to remove bad ids

    my $cmd = "perl $bin/filter_sam_from_prinseq.pl --sam_file=$bam_file --bad_list=$output_dir/$name\_bad_ids.out --out_file=$output_dir/$name\_filtered.sam";
    $self->_run_cmd($cmd);

    my $cmd = "$samtools view -S -b $output_dir/$name\_filtered.sam > $output_dir/$name\_filtered.bam";
    $self->_run_cmd($cmd);


    my $count = $self->_run_cmd("$samtools view $output_dir/$name\_filtered.bam | cut -f1 | uniq | wc -l");
    chomp $count;

    # Blitz the sam file
    $self->_run_cmd("rm $output_dir/$name\_filtered.sam");

#    my $cmd = "perl $bin/sam2fasta.pl --input=$output_dir/$name\_filtered.bam --combine_mates=0 --paired=1 --output_file=$output_dir/$name.fastq";
#    $self->_run_cmd($cmd);
#    `rm $output_dir/$name\_filtered.sam`; 
    return {
        count => $count,
        file => "$output_dir/$name\_filtered.bam"
    }

}

=head2 run_cmd

 Title   : sam2Fasta
 Usage   : my $fastas = $LGTSeek->sam2Fasta({'input' => '/path/to/file.bam'...})
 Function: Convert a bam/sam file to a fasta file
 Returns : a list of fasta/fastq files
 Args    : sam or bam file to convert

=cut
sub sam2Fasta {
    my($self, $config) = @_;
    my $bin = $self->{ergatis_bin};

    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    my $outfile;
    my($name,$path,$suff) = fileparse($config->{input},(qr/.bam$||.sam$||.sam.gz$/));
    my $cmd = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$config->{input}";
    if($config->{fastq}) {
        $outfile = "$output_dir/$name.fastq";
        $cmd .= " --fastq=1 --output_file=$outfile";
    }
    else {
        $outfile = "$output_dir/$name.fasta";
        $cmd .= " --fastq=0 --output_file=$outfile";
    }
    if($config->{combine_mates}) {
        $cmd .= " --combine_mates=0";
    }
    if($config->{paired}) {
        $cmd .= " --paired=1";
    }
    
    $self->_run_cmd($cmd);

    return $outfile;
}

=head2 run_cmd

 Title   : _run_cmd
 Usage   : *PRIVATE*
 Function: Run a unix command and fail if something goes wrong
 Returns : void
 Args    : Command to run

=cut
sub _run_cmd {

    my($self, $cmd) = @_;

    my $res = `$cmd`;
    if($?) {
        print STDERR "$cmd\n\n$?";
        die "$cmd died with message:\n$res\n\n";
    }
    return $res;
   # print "$cmd\n";
}

=head2 downloadSRA

 Title   : downloadSRA
 Usage   : $lgtseek->downloadSRA(({'experiment_id'} = 'SRX01234'})
 Function: Download sra files from the sequence read archive
 Returns : A list of the downloaded file paths
 Args    : 

=cut
sub downloadSRA {
    my ($self, $config) = @_;

    # Check for all the aspera related options
    $self->{aspera_host} = $config->{aspera_host} ? $config->{aspera_host} : $self->{aspera_host};
    $self->{aspera_user} = $config->{aspera_user} ? $config->{aspera_user} : $self->{aspera_user};    
    $self->{aspera_params} = $config->{aspera_params} ? $config->{aspera_params} : $self->{aspera_params};
    $self->{aspera_path} = $config->{aspera_path} ? $config->{aspera_path} : $self->{aspera_path};
    $self->{aspera_rate} = $config->{aspera_rate} ? $config->{aspera_rate} : $self->{aspera_rate};
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->{aspera_path} = $self->{aspera_path} ? $self->{aspera_path} : '~/.aspera/connect/';

    my $retry_attempts = $config->{retry_attempts} ? $config->{retry_attempts} : 10;
    if(!$self->{aspera_host}) {
        $self->{aspera_host} = 'ftp-private.ncbi.nlm.nih.gov';
    }
    if(!$self->{aspera_user}) {
        $self->{aspera_user} = 'anonftp';
    }
    if(!$self->{aspera_rate}) {
        $self->{aspera_rate} = '200M';
    }

    if(!$self->{aspera_path}) {
        die "Need to specify an aspera_path (where is aspera installed) in order to download from the sra\n";
    }
    if(!$self->{output_dir}) {
        die "Need to specify an output_dir in order to download from the sra\n";
    }
    
    my $prefix;
    my $thousand;
    my $output_dir = $self->{output_dir};
    my $path_to_file;


    # We can pass an experiment_id, run_id or a full path to download
    if($config->{experiment_id}) {
        my $exp_id = $config->{experiment_id};
        $exp_id =~ /^((.{3}).{3}).*/;
        $prefix = $2;
        $thousand = $1;
        $path_to_file = "/sra/sra-instant/reads/ByExp/litesra/$prefix/$thousand/$exp_id";
        $output_dir = "$output_dir/$prefix/$thousand/";
    }
    if($config->{run_id}) {
        $config->{run_id} =~ /^((.{3}).{3}).*/;
        $prefix = $2;
        $thousand = $1;
        $path_to_file = "/sra/sra-instant/reads/ByRun/litesra/$prefix/$thousand/$config->{run_id}";
        $output_dir = "$output_dir/$prefix/$thousand/";
    }
    elsif($config->{path}) {
        $path_to_file = $config->{path};
    }

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");


    my $cmd_string = "$self->{aspera_path}/bin/ascp -QTd -l$self->{aspera_rate} -i $self->{aspera_path}/etc/asperaweb_id_dsa.putty $self->{aspera_user}\@$self->{aspera_host}:$path_to_file $output_dir -L $output_dir -o Overwrite=diff 2>&1";


    #Retry the download several times just incase.
    my $retry = 1;

    my $retries = 0;
    while($retry) {
        
        # Doing this echo y to ensure we accept any certs.
        my $out = $self->_run_cmd("echo y | $cmd_string");

        # We can actually exit non-0 and still succeed if the 
        if($out =~ /Error/)  {
            print STDERR "Had a problem downloading $self->{aspera_host}:$path_to_file to $output_dir\n";
            print STDERR "$cmd_string";
            if($retries < $retry_attempts) {
                $retries++;
                sleep $retries*2; # Sleep for 2 seconds per retry.
            }
            else {
                print STDERR "Retries exhausted. Tried $retries times.\n$cmd_string";
                exit(1);
            }
        }
        else {
            $retry = 0;
            print STDERR "$? $out Successfully downloaded $self->{aspera_host}:$path_to_file to $output_dir\n";
        }
     }   

    my @files = `find $output_dir -name *.sra`;
    return \@files;
}

=head2 downloadCGHub

 Title   : downloadCGHub
 Usage   : $lgtseek->downloadCGHub(({'analysis_id'} = '00007994-abeb-4b16-a6ad-7230300a29e9'})
 Function: Download TCGA bam files from the CGHub
 Returns : A list of the downloaded file paths
 Args    : 

=cut

sub downloadCGHub {
    my ($self, $config) = @_;

    # Check for all the genetorrent related options
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->{genetorrent_path} = $self->{genetorrent_path} ? $self->{genetorrent_path} : '/opt/genetorrent/bin';
    $self->{cghub_key} = $config->{cghub_key} ? $config->{cghub_key} : $self->{cghub_key};

    my $retry_attempts = $config->{retry_attempts} ? $config->{retry_attempts} : 10;
    if(!$self->{cghub_key}) {
#        $self->{aspera_user} = 'anonftp';
         die "Need to specify the path to the cghub_key\n";
    }

    if(!$self->{genetorrent_path}) {
        die "Need to specify an genetorrent_path (where is genetorrent/gtdownload installed) in order to download from CGHub\n";
    }
    if(!$self->{output_dir}) {
        die "Need to specify an output_dir in order to download from CGHub\n";
    }

    my $output_dir = $self->{output_dir};

    # We can pass an analysis_id (UUID), URI, .xml .gto to download. They can all be called analysis_id and will work the same way.
    my $download = $config->{analysis_id};

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");

    my $cmd_string = "$self->{genetorrent_path}/gtdownload -d $download -p $output_dir -c $self->{cghub_key}";
    $self->_run_cmd($cmd_string);

    #Retry the download several times just incase.
#    my $retry = 1;

#    my $retries = 0;
#    while($retry) {

        # Doing this echo y to ensure we accept any certs.
#        my $out = $self->_run_cmd("echo y | $cmd_string");

        # We can actually exit non-0 and still succeed if the 
#        if($out =~ /Error/)  {
#            print STDERR "Had a problem downloading $download to $output_dir\n";
#            print STDERR "$cmd_string";
#            if($retries < $retry_attempts) {
#                $retries++;
#                sleep $retries*2; # Sleep for 2 seconds per retry.
#            }
#            else {
#                print STDERR "Retries exhausted. Tried $retries times.\n$cmd_string";
#                exit(1);
#            }
#        }
#        else {
#            $retry = 0;
#            print STDERR "$? $out Successfully downloaded $download to $output_dir\n";
#        }
#     }

    my @bam = `find $output_dir -name *.bam`;
    my @bai = `find $output_dir -name *.bai`;
    my @gto = `find $output_dir -name *.gto`;

    my $files = {
        'bam_files' => \@bam,
        'bai_files' => \@bai,
        'gto_files' => \@gto
    };
    
    return $files;
}


=head2 dumpFastq

 Title   : dumpFastq
 Usage   : $lgtseek->dumpFastq(({'sra_file' => 'SRX01234'})
 Function: Run the sratoolkit program dump-fastq
 Returns : An object with a path and basename of the output as well as a list of output files
 Args    : An object with element 'sra_file' and optionally the path to the sratoolkit install

=cut

sub dumpFastq {
    my($self,$config) = @_;

    $self->{sratoolkit_path} = $config->{sratoolkit_path} ? $config->{sratoolkit_path} : $self->{sratoolkit_path};

    # If we don't have a path provided we'll hope it's in our path.
    my $fastqdump_bin = $self->{sratoolkit_path} ? "$self->{sratoolkit_path}/fastq-dump" : "fastq-dump";

    $config->{sra_file} =~ s/\/\//\//g;
    # Need to pull the version of the sratoolkit to determine if we need the --split-3 parameter.
    my $ret = `$fastqdump_bin -V`;
    my $version;
    my $cutoff_version;
    if($ret =~ /fastq-dump : ([\d.]+)/) {
        $version = version->parse($1);
        $cutoff_version = version->parse('2.1.0');
    }
    else {
        die "$? $ret $fastqdump_bin\n";
    }
    if($version > $cutoff_version) {

        $fastqdump_bin .= " --split-3 ";
    }
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if(!$self->{output_dir}) {
        die "Need to specify an output_dir in order to download from the sra\n";
    }

    my ($name,$path,$suff) = fileparse($config->{sra_file},'.sra');
    chomp $name;
    my $res = $self->_run_cmd("find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
    my @files = split(/\n/,$res);

    if(!@files && !$config->{overwrite}) {
        my $cmd = "$fastqdump_bin -O $self->{output_dir} $config->{sra_file}";
        $self->_run_cmd($cmd);
        my $res = $self->_run_cmd("find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
        @files = split(/\n/,$res);
    }

    print STDERR "@files\n";

    my $retval = {
        'files' => \@files,
        'path' => $self->{output_dir},
        'basename' => $name,
        'paired_end' => 0
    };
    if($files[0] =~ /_\d.fastq/) {
        $retval->{paired_end} = 1;
    }
    return $retval;
}

=head2 runBWA

 Title   : runBWA
 Usage   : $lgtseek->runBWA(({'base' => 'SRX01234','path' => '/path/to/files'})
 Function: Run bwa using the lgt_bwa wrapper
 Returns : The path to the bam file
 Args    : The input fasq/bam files and references which can be done a few different ways:

           # For files like /path/to/files/SRR01234_1.fastq and /path/to/files/SRR01234_2.fastq
           {'input_dir' => '/path/to/files/',
            'input_base' => 'SRR01234',
            'reference' => '/path/to/references/hg19.fa'
           }
         
           # For bam files and a list of references
           {'input_bam' => '/path/to/files/SRR01234.bam',
            'reference_list' => '/path/to/references/all_refs.list'
           }

=cut

sub runBWA {
    my($self,$config) = @_;

    $self->{ergatis_bin} = $config->{ergatis_bin} ? $config->{ergatis_bin} : $self->{ergatis_bin};
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    # Check for a bwa path. If we don't have one we'll just hope it's in our global path.
    $self->{bwa_path} = $config->{bwa_path} ? $config->{bwa_path} : $self->{bwa_path};
    $self->{bwa_path} = $self->{bwa_path} ? $self->{bwa_path} : 'bwa';
    
    $self->_run_cmd("mkdir -p $output_dir");

    my $conf = {
        num_aligns => 0,
        bwa_path => $self->{bwa_path},
        output_dir => $output_dir
    };


    # Build the command string;
    my @cmd = ("$self->{ergatis_dir}/lgt_bwa --num_aligns=0 --bwa_path=$self->{bwa_path}");

    my $suff = '.sam';
    if($config->{output_bam}) {
        $suff = '.bam';
        $conf->{output_bam} = 1;
#        push(@cmd, '--output_bam=1');
    }

    my $basename = $config->{input_base};

    # Handle making up the lgt_bwa command with a bam file
    if($config->{input_bam}) {
        my ($name,$path,$suff) = fileparse($config->{input_bam},'.bam');
        $basename = $name;
        $conf->{input_base}=$basename;
        $conf->{input_bam} = $config->{input_bam};
#        push(@cmd,"--input_bam=$config->{input_bam}");
    }
    elsif($config->{input_dir} && $config->{input_base}) {
        $conf->{input_dir} = $config->{input_dir};
        $conf->{input_base} = $config->{input_base};
#        push(@cmd,"--input_dir=$config->{input_dir} --input_base=$config->{input_base}");
    }
    else {
        die "Must provide either a value to either input_bam or to input_base and input_dir\n";
    }

    my $pre = '';
    if($config->{reference}) {
        my ($name,$dir,$suff) = fileparse($config->{reference},qr/\.[^\.]+/);
        $pre = "$name\_";
        $conf->{ref_file} = $config->{reference};
#       push(@cmd,"--ref_file=$config->{reference");
    }
    elsif($config->{reference_list}) {
        $conf->{ref_file_list} = $config->{reference_list};
#        push(@cmd,"--ref_file_list=$config->{reference_list}");  
    }
    else {
        die "Must provide a value for either reference or reference_list to run bwa\n";
    }

    $conf->{overwrite} = $config->{overwrite};
#    push(@cmd, "--output_dir=$output_dir");
    map {
        $conf->{$_} = $config->{other_opts}->{$_};
    } keys %{$config->{other_opts}};
#    push(@cmd, $config->{other_opts});
    $conf->{run_lca} = $config->{run_lca};
    $conf->{lgtseek} = $self;
    $conf->{cleanup_sai} = $config->{cleanup_sai};
    $conf->{out_file} = $config->{out_file};
    print STDERR "about to call runBWA\n";
    LGTbwa::runBWA($conf);

    # Maybe should check if this is valid.
#    my $res = $self->_run_cmd($cmd_string);
    if($config->{run_lca}) {
    
    }
    else {
    
        my @files = split(/\n/,$self->_run_cmd("find $output_dir -name '*$pre$basename$suff'"));
        map {chomp $_;} @files; 
        print STDERR join(" ",@files);
        print STDERR "\n";
        return \@files;
    }
}

=head2 bwaPostProcess

 Title   : bwaPostProcess
 Usage   : $lgtseek->bwaPostProcess(({'donor_bams' => \@donors,'host_bams' => \@hosts})
 Function: Classify the results of a short read mapping (UM, UU, MM, UM_UM etc.)
 Returns : An object with counts of the different classes as well as the path to bam files 
           containing these reads.
 Args    : An object with donor and optionally host bam files.

=cut

sub bwaPostProcess {
    my ($self, $config) = @_;

    my $retval;
    # Do we have both donor and host bams?
    if($config->{donor_bams} && $config->{host_bams}) {
        $retval = $self->_bwaPostProcessDonorHostPaired($config);
    }
    
    return $retval;
}

=head2 _bwaPostProcessDonorPaired

 Title   : _bwaPostProcessDonorPaired
 Function: Classify the results of a short read mapping (UM, UU, MM, etc.)
 Returns : An object with counts of the different classes as well as the path to bam files 
           containing these reads.
 Args    : An object with donor and optionally host bam files.

=cut
sub _bwaPostProcessDonorPaired {
    my ($self, $config) = @_;

    my @donor_fh;
    my @donor_head;

    $self->{samtools_bin} = $self->{samtools_bin} ? $self->{samtools_bin} : 'samtools';
    my $samtools = $self->{samtools_bin};
    my $prefix = $config->{output_prefix} ? "$config->{output_prefix}_" : '';

    # Open all the donor files
    map {
        print STDERR "Opening $_\n";

        if($_ =~ /.bam$/) {
            push(@donor_head, `$samtools view -H $_`);
            open(my $fh, "-|", "$samtools view $_");
            push(@donor_fh,$fh);
        }
        elsif($_ =~ /.sam.gz$/) {
            push(@donor_head, `zcat $_ | $samtools view -H -S -`);
            open(my $fh, "-|", "zcat $_ | $samtools view -S -");
            push(@donor_fh,$fh);
        }
    } @{$config->{donor_bams}};

    my $class_to_file_name = {
        'lgt_donor' => "$self->{output_dir}/".$prefix."lgt_donor.bam",
        'microbiome_donor' => "$self->{output_dir}/".$prefix."microbiome.bam"
    };

    
    # Check if these files exist already. If they do we'll skip regenerating them.
    my $files_exist = 1;

    map {
        if(! -e $class_to_file_name->{$_}) {
            $files_exist = 0;
        }
    } keys %$class_to_file_name;

    my $class_counts = {
        'lgt' => undef,
        'microbiome' => undef
    };

    # If the files are already there and we aren't being forced to overwrite, we'll
    # just get the counts and return
    if($files_exist && !$config->{overwrite}) {
        
        map {
            $_ =~ /^*(.*)\_\w+$/;
            my $class = $1;
            if(!$class_counts->{$class}) {
                my $count = $self->_run_cmd("$samtools view $class_to_file_name->{$_} | wc -l");
                chomp $count;
                $class_counts->{$class} = $count;
                print STDERR "$count for $class\n";
            }
        } keys %$class_to_file_name;
        return {
            counts => $class_counts,
            files => $class_to_file_name
        };
    }
    # Here are a bunch of file handles we'll use later.
    open(my $lgtd,  "| $samtools view -S -b -o $self->{output_dir}/".$prefix."lgt_donor.bam -") or die "Unable to open\n";
    open(my $microbiome_donor,"| $samtools view -S -b -o $self->{output_dir}/".$prefix."microbiome.bam -") or die "Unable to open\n";

    my $class_to_file = {

    };

    my $more_lines = 1;

    while($more_lines) {
        
        my @donor_lines;
        
        my $dr1_line;
        my $dr2_line;

        my $obj = $self->_getPairedClass({fhs => \@donor_fh});
        my $class = $obj->{class};
        $more_lines = $obj->{more_lines};
        $dr1_line = $obj->{r1_line};
        $dr2_line = $obj->{r2_line};


#        if($class_to_file->{$classes_both->{$paired_class}."_donor"}) {
#            print STDERR "Printing to ".$classes_both->{$paired_class}."_donor $paired_class\n$dr1_line$dr2_line";
#            print {$class_to_file->{$classes_both->{$paired_class}."_donor"}} "$dr1_line\n$dr2_line\n";
#        }
#        if($class_to_file->{$classes_both->{$paired_class}."_host"}) {
#            print {$class_to_file->{$classes_both->{$paired_class}."_host"}} "$hr1_line\n$hr2_line\n";
#        }

#        if($classes_both->{$paired_class} eq 'lgt') {
#            print STDERR "Processing $hr1_line$hr2_line$dr1_line$dr2_line";
#        }
#        if($classes_both->{$paired_class}) {
#            $class_counts->{$classes_both->{$paired_class}}++;
#            my @lines;
#            map {push(@lines, "$_: $class_counts->{$_}")} keys %$class_counts;
#            map {print  "\r$_: $class_counts->{$_}"} keys %$class_counts;
#                   print "\r".join(' ',@lines);
                
                #print STDERR "Line $line_num donor class: $dclass\nHost class: $hclass\nCombined class: $paired_class\n";
 #       }
    }
}

sub _getPairedClass {
    my ($self, $config) = @_;

    my $fhs = $config->{fhs};
    my $more_lines = 1;
    
    my $r1_class;
    my $r1_line;
    my $r2_class;
    my $r2_line;
    # Next establish the class of the donor read
    foreach my $fh (@$fhs) {
        my $r1 = <$fh>;
        my $r2 = <$fh>;
        if($config->{strip_xa}) {
            chomp $r1;
            $r1 =~ s/\tXA:Z\S+$//;
            chomp $r2;
            $r2 =~ s/\tXA:Z\S+$//;
        }  
#                print STDERR "Processing $hr1$hr2$dr1$dr2";
        # Should check if these ended at the same time?
        if(!$r1 || !$r2) {
            $more_lines = 0;
            last;
        }
        
        my $r1_flag = $self->_parseFlag((split(/\t/,$r1))[1]);
#            my $dr2_flag = $self->_parseFlag((split(/\t/,$dr2))[1]);
        if(!$r1_flag->{'qunmapped'}) {
            $r1_line = $r1;
            $r1_class = 'M';
        }
        elsif(!$r1_class) {
            $r1_line = $r1;
            $r1_class = 'U';
        }
        if(!$r1_flag->{'munmapped'}) {
            $r2_line = $r2;
            $r2_class = 'M';
        }            
        elsif(!$r2_class) {
            $r2_line = $r2;
            $r2_class = 'U';
        }
    }

    my $class = "$r1_class$r2_class";

    return {
        class => $class,
        r1_line => $r1_line,
        r2_line => $r2_line,
        more_lines => $more_lines
    }
}

sub _bwaPostProcessDonorHostPaired {
    my ($self, $config) = @_;

    $self->{samtools_bin} = $self->{samtools_bin} ? $self->{samtools_bin} : 'samtools';
    my $samtools = $self->{samtools_bin};
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    my $classes_each = {
        'MM' => 'paired',
        'UM' => 'single',
        'MU' => 'single',
        'UU' => 'none'
    };

    # Classes have the donor on the left and the host on the right
    my $classes_both = {
        'UM_MU' => 'lgt',
        'MU_UM' => 'lgt',
        'UM_MM' => 'integration_site_host',
        'MU_MM' => 'integration_site_host',
        'UU_UM' => 'integration_site_host',
        'UU_MU' => 'integration_site_host',
        'MM_UM' => 'integration_site_donor',
        'MM_MU' => 'integration_site_donor',
        'MU_UU' => 'integration_site_donor',
        'UM_UU' => 'integration_site_donor',
        'MM_UU' => 'microbiome',
        'UU_MM' => 'host',
        'UU_UU' => 'no_map',
        'MM_MM' => 'all_map',
        'UM_UM' => 'single_map',
        'MU_MU' => 'single_map'
    };

    my $prefix = $config->{output_prefix} ? "$config->{output_prefix}_" : '';

    my $class_to_file_name = {
        'lgt_donor' => "$self->{output_dir}/".$prefix."lgt_donor.bam",
        'lgt_host' => "$self->{output_dir}/".$prefix."lgt_host.bam",
        'integration_site_donor_donor' => "$self->{output_dir}/".$prefix."integration_site_donor_donor.bam",
        'integration_site_donor_host' => "$self->{output_dir}/".$prefix."integration_site_donor_host.bam",
        'microbiome_donor' => "$self->{output_dir}/".$prefix."microbiome.bam",
    };

    
    # Check if these files exist already. If they do we'll skip regenerating them.
    my $files_exist = 1;

    map {
        if(! -e $class_to_file_name->{$_}) {
            $files_exist = 0;
        }
    } keys %$class_to_file_name;

    my $class_counts = {
        'lgt' => undef,
        'integration_site_host' => undef,
        'integration_site_donor' => undef,
        'microbiome' => undef,
        'host' => undef,
        'no_map' => undef,
        'all_map' => undef,
        'single_map' => undef
    };

    # If the files are already there and we aren't being forced to overwrite, we'll
    # just get the counts and return
    if($files_exist && !$config->{overwrite}) {
        
        map {
            $_ =~ /^*(.*)\_\w+$/;
            my $class = $1;
            if(!$class_counts->{$class}) {
                my $count = $self->_run_cmd("$samtools view $class_to_file_name->{$_} | cut -f1 | uniq | wc -l");
                chomp $count;
                $class_counts->{$class} = $count*1;
                print STDERR "$count for $class\n";
            }
        } keys %$class_to_file_name;
        my $total = $self->_run_cmd("$samtools view ".$config->{donor_bams}[0]." | cut -f1 | uniq | wc -l");
        chomp $total;
        $class_counts->{total} = $total*1;
        return {
            counts => $class_counts,
            files => $class_to_file_name
        };
    }
    # Here are a bunch of file handles we'll use later.
    print STDERR "$self->{output_dir}/".$prefix."lgt_donor.bam\n";
    open(my $lgtd,  "| $samtools view -S -b -o $self->{output_dir}/".$prefix."lgt_donor.bam -") or die "Unable to open\n";
    open(my $lgth, "| $samtools view -S -b -o $self->{output_dir}/".$prefix."lgt_host.bam -") or die "Unable to open\n";
    open(my $int_site_donor_d, "| $samtools view -S -b -o $self->{output_dir}/".$prefix."integration_site_donor_donor.bam -") or die "Unable to open\n";
    open(my $int_site_donor_h, "| $samtools view -S -b -o $self->{output_dir}/".$prefix."integration_site_donor_host.bam -") or die "Unable to open\n";
    open(my $microbiome_donor,"| $samtools view -S -b -o $self->{output_dir}/".$prefix."microbiome.bam -") or die "Unable to open\n";


    my $class_to_file = {
        'lgt_donor' => $lgtd,
        'lgt_host' => $lgth,
        'integration_site_donor_donor' => $int_site_donor_d,
        'integration_site_donor_host' => $int_site_donor_h,
        'microbiome_donor' => $microbiome_donor
    };

    my @donor_fh;
    my @host_fh;
    my @donor_head;
    my @host_head;

    # Open all the donor files
    map {
        print STDERR "Opening $_\n";

        if($_ =~ /.bam$/) {
            push(@donor_head, `$samtools view -H $_`);
            open(my $fh, "-|", "$samtools view $_");
            push(@donor_fh,$fh);
        }
        elsif($_ =~ /.sam.gz$/) {
            push(@donor_head, `zcat $_ | $samtools view -H -S -`);
            open(my $fh, "-|", "zcat $_ | $samtools view -S -");
            push(@donor_fh,$fh);
        }
    } @{$config->{donor_bams}};
    
    # Open all the host files
    map {
        print STDERR "Opening $_\n";
        if($_ =~ /.bam$/) {
            push(@host_head, `$samtools view -H $_`);
            open(my $fh, "-|", "$samtools view $_");
            push(@host_fh,$fh);
        }
        elsif($_ =~ /.sam.gz$/) {
            push(@host_head, `zcat $_ | $samtools view -H -S -`);
            open(my $fh, "-|", "zcat $_ | $samtools view -S -");
            push(@host_fh,$fh);
        }
    } @{$config->{host_bams}};


    # Prime the files with headers.
    map {
        if($_ =~ /_donor$/) {
            print "Printing header to $_ donor file\n";
            print { $class_to_file->{$_}} join('',@donor_head);
#            close $class_to_file->{$_};
        }
        elsif($_ =~ /_host$/) {
            print "Printing header to $_ host file\n";
            print {$class_to_file->{$_}} join('',@host_head);
#            close $class_to_file->{$_};
        }
    }keys %$class_to_file;
 #   exit;
    my $more_lines = 1;

    my $line_num =0;

    while($more_lines) {
        
        my @donor_lines;
        my @host_lines;
        
        # Get the class of the host mappings
        my $obj = $self->_getPairedClass({fhs => \@host_fh});
        my $hclass = $obj->{class};
        $more_lines = $obj->{more_lines};
        my $hr1_line = $obj->{r1_line};
        my $hr2_line = $obj->{r2_line};

        # Get the class of the donor mappings
        my $obj = $self->_getPairedClass({fhs => \@donor_fh});
        my $dclass = $obj->{class};
        $more_lines = $obj->{more_lines};
        my $dr1_line = $obj->{r1_line};
        my $dr2_line = $obj->{r2_line};

        if($more_lines) {
            my $paired_class = "$dclass\_$hclass";

            # print the donor lines to the donor file (if we are keeping this output file)
            if($class_to_file->{$classes_both->{$paired_class}."_donor"}) {
                print {$class_to_file->{$classes_both->{$paired_class}."_donor"}} "$dr1_line\n$dr2_line\n";
            }
            
            # print the host lines to the host file (if we are keeping this output file)
            if($class_to_file->{$classes_both->{$paired_class}."_host"}) {
                print {$class_to_file->{$classes_both->{$paired_class}."_host"}} "$hr1_line\n$hr2_line\n";
            }
            
            # Increment the count for this class
            if($classes_both->{$paired_class}) {
                $class_counts->{$classes_both->{$paired_class}}++;
            }
            # Increment the total count
            $line_num ++;
        }
    }

    # Close up the file handles
    map {
        if($_ =~ /_donor$/) {
            print STDERR "closing $_ donor file\n";
            close $class_to_file->{$_};
        }
        elsif($_ =~ /_host$/) {
            print STDERR "closing $_ host file\n";
            close $class_to_file->{$_};
        }
    }keys %$class_to_file;

    # Set the total
    $class_counts->{total} = $line_num;

    # Return the files and the counts
    return {
        counts => $class_counts,
        files => $class_to_file_name
    };
}

sub blast2lca {
    my ($self,$config) = @_;

    my $gi2tax = $self->getGiTaxon({});

    open IN, "<$config->{blast}" or die "Unable to open $config->{blast}\n";
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    
    my ($name,$directories,$suffix) = fileparse($config->{blast},qr/\.[^.]*/);
    open OUT, ">$output_dir/$name\_lca.out" or die "Unable to open $output_dir/$name\_lca.out";
    print OUT join("\t", ('read_id','best_evalue','lca','paired_lca','liberal_lca'));
    print OUT "\n";
    my $hits_by_readname = {};
    my $id;
    my $evalue;
    my $lca;

    my @vals;
    while (<IN>) {
        my @fields = split(/\t/);
        my $new_id = $fields[0];
        if($fields[0] =~ /^(.*)\/\d/) {
            $new_id = $1;
        }
        if($new_id eq $id && $fields[10] > $evalue) {
            next;
        }
        my $taxon =  $gi2tax->getTaxon($fields[1]);
        $taxon->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;


#        if(&filter_genera($genera)) {
#            next;
#        }
        if(!defined($id)) {
            $id = $new_id;
            $evalue = $fields[10];
            $lca = $taxon->{lineage};
            $hits_by_readname->{$fields[0]} = {$taxon->{taxon_id}} => $taxon->{lineage};
        }elsif($new_id eq $id && $fields[10] == $evalue) {
            $hits_by_readname->{$fields[0]}->{$taxon->{taxon_id}} = $taxon->{lineage};
            my $newlca = &_find_lca([$lca, $taxon->{lineage}]);
            $lca = $newlca;
        }
        elsif($new_id ne $id) {

            @vals = ($id,$evalue,$lca);
            my $newvals = &_prep_lca_print($hits_by_readname,\@vals);
            

            print OUT join("\t",@$newvals);
            print OUT "\n";
            $id = $new_id;
            $evalue = $fields[10];
            $lca = $taxon->{lineage};
        }
    }
    # Print out last line

    @vals = ($id,$evalue,$lca);
    my $newvals = &_prep_lca_print($hits_by_readname,\@vals);
    print OUT join("\t",@$newvals);
    print OUT "\n";
    return "$output_dir/$name\_lca.out";
}

sub _prep_lca_print {
    my $hits_by_readname = shift;
    my $vals = shift;
    # Find the LCA of only mated hits
    my @reads = keys %$hits_by_readname;
    if(scalar @reads > 1) {
        
        
        # First we'll find the LCA requiring both reads to match the same
        # taxon id.
        my @good_lineages;
        my @read1_lineages;
        map {
            if($hits_by_readname->{$reads[1]}->{$_}) {
                push(@good_lineages,$hits_by_readname->{$reads[1]}->{$_});
            }
            push(@read1_lineages,$hits_by_readname->{$reads[0]}->{$_});
        }keys %{$hits_by_readname->{$reads[0]}};
        
        my @read2_lineages;
        map {
            push(@read2_lineages,$hits_by_readname->{$reads[1]}->{$_});
        }keys %{$hits_by_readname->{$reads[1]}};
        my $paired_lca = &_find_lca(\@good_lineages);
        push(@$vals,$paired_lca);
        
        # Next we'll find the LCA where if either LCA is a substring of the
        # other we'll take the longer one (more specific).
        my $read1_lca = &_find_lca(\@read1_lineages);
        my $read2_lca = &_find_lca(\@read2_lineages);
        
        my $liberal_lca = &_find_lca([$read1_lca,$read2_lca]);

        if($read1_lca =~ $read2_lca) {
            $liberal_lca = $read2_lca;
        }
        elsif($read2_lca =~ $read1_lca) {
            $liberal_lca = $read1_lca;
        }

        push(@$vals,$liberal_lca);
        $hits_by_readname = {};
    }
    else {
        push(@$vals,('',''));
    }
    return $vals;
}


sub _find_lca {
    my $lineages = shift;

    # prime it
    my @lca = split(';', $lineages->[0]);

    foreach my $l (@$lineages) {
        my $newlca = [];
        my @lineage = split(';',$l);
        for( my $i = 0; $i < @lineage;$i++) {
            if($lca[$i] eq $lineage[$i]) {
                push(@$newlca, $lineage[$i]);
            }
            else {
                last;
            }   
        }
        @lca = @$newlca;
    }
    if(! scalar @lca) {
#        print STDERR "Had no lca @$lineages\n";
    }
    #print STDERR join(";",@lca);
    #print STDERR "\n";
    return join(';',@lca);
}

sub bestBlast2 {
    my ($self,$config) = @_;
    my $fasta;
    # Convert bams to fasta
    if($config->{'bam'}) {
        my $newfasta = $self->sam2Fasta({
            input => $config->{'bam'},
            output_dir => $config->{output_dir},
            paired => 1});
        $fasta = $newfasta;
    }
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    if($config->{'fasta'}) {
        $fasta = $config->{fasta};
    }
#    my $blast_bin = '';
#    if($config->{blast_bin}) {
        
        #$blast_bin = $config->{"blast_bin"}." -d $config->{'blast_db'} -e 1e-5 -m8 -T F ";
#    }  

    my ($name,$directories,$suffix) = fileparse($fasta,qr/\.[^.]*/);
    open OUT, ">$self->{output_dir}/$name\_filtered_blast.list" or die "Unable to open $self->{output_dir}/$name\_filtered_blast.list\n";  
    my @outputs;
#    foreach my $file (@fastas) {
        my $files = LGTBestBlast::filterBlast({
            blast_bin => $config->{blast_bin},
            fasta => $fasta,
            db => $config->{db},
            gitaxon => $self->getGiTaxon({}),
            lineage1 => $config->{lineage1},
            lineage2 => $config->{lineage2},
            output_dir => "$self->{output_dir}"});
        map {
            print OUT "$files->{$_}\n";
        } keys %$files;

#        push(@outputs,$files);
 #   }
    close OUT;

    $files->{list_file} = "$self->{output_dir}/$name\_filtered_blast.list";
    return $files;
}
=head2 bestBlast

 Title   : bestBlast
 Usage   : $lgtseek->bestBlast(({'inputs' => \@bam_files, 'ref' => '/path/to/ref/db'})
 Function: Run blast (or megablast) against a reference database to find best blast hits. Then
           determine if mates look like valid LGT's
 Returns : A file with all the LGT's hit information
 Args    : An object with the input bame files and the reference database.

=cut
sub bestBlast {
    my ($self,$config) = @_;
    
    my @fastas;
    # Convert bams to fasta
    if($config->{'bams'}) {
    
    }    
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    if($config->{'fastas'}) {
        push(@fastas,@{$config->{'fastas'}});
    }
    my $blast_bin = '';
    if($config->{blast_bin}) {
        $blast_bin = $config->{"blast_bin"}." -d $config->{'blast_db'} -e 1 -m8 -T F ";
    }        
    my @overall_files;
    my @lineage1_files;
    my @lineage2_files;
    foreach my $fasta (@fastas) {
        
        my $cmd = "perl -I $self->{bin_dir} $self->{bin_dir}/filter_lgt_best_hit.pl".
        " --input=$fasta".
        " --blast_bin=\'$blast_bin\'".
        " --taxon_dir=".$self->{taxon_dir}.
        " --trace_mapping=".$config->{trace_mapping}.
        " --lineage1=".$config->{lineage1}.
        " --lineage2=".$config->{lineage2}.
        " --filter_lineage=".$config->{filter_lineage}.
        " --dbhost=".$self->{taxon_host}.
        " --idx_dir=".$self->{taxon_idx_dir}.
        " --output_dir=".$self->{output_dir};
        
        
        $self->_run_cmd($cmd);
        my ($name,$directories,$suffix) = fileparse($fasta,qr/\.[^.]*/);
        
        push(@overall_files,"$self->{output_dir}/$name\_overall.out");
        push(@lineage1_files,"$self->{output_dir}/$name\_lineage1.out");
        push(@lineage2_files,"$self->{output_dir}/$name\_lineage2.out");
    }
    
    # Should merge the three file types here.
    my $cmd = "cat ".join(' ',@overall_files)." > $self->{output_dir}/all_overall.out";
    $self->_run_cmd($cmd);
    my $cmd = "cat ".join(' ',@lineage1_files)." > $self->{output_dir}/all_lineage1.out";
    $self->_run_cmd($cmd);
    my $cmd = "cat ".join(' ',@lineage2_files)." > $self->{output_dir}/all_lineage2.out";
    $self->_run_cmd($cmd);

    open OUT, ">$self->{output_dir}/filtered_blast.list" or die "Unable to open $self->{output_dir}/filtered_blast.list\n";
    print OUT "$self->{output_dir}/all_overall.out\n$self->{output_dir}/all_lineage1.out\n$self->{output_dir}/all_lineage2.out";
    close OUT;

    return {
        overall_blast => "$self->{output_dir}/all_overall.out",
        lineage1_blast => "$self->{output_dir}/all_lineage1.out",
        lineage2_blast => "$self->{output_dir}/all_lineage2.out",
        list_file => "$self->{output_dir}/filtered_blast.list"
    };
}

=head2 bestBlast

 Title   : bestBlast
 Usage   : $lgtseek->runLgtFinder(({'inputs' => \@files})
 Function: Run blast (or megablast) against a reference database to find best blast hits. Then
           determine if mates look like valid LGT's
 Returns : A file with all the LGT's hit information
 Args    : An object with the input bame files and the reference database.

=cut
sub runLgtFinder {
    my ($self,$config) = @_;

    LGTFinder::findLGT($config);

#    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

#    my $cmd = "perl $self->{ergatis_dir}/lgt_finder_dr.pl".
#        " --input_file_list=$config->{input_file_list}".
#        " --output_prefix=$config->{output_prefix}".
#        " --output_dir=$self->{output_dir}".
#        " --ref_lineage=$config->{ref_lineage}";
        
#    $self->_run_cmd($cmd);

    my $pref = $config->{output_prefix} ? $config->{output_prefix} : 'lgt_finder';

    my $valid_count = $self->_run_cmd("grep ';$config->{lineage1};' $self->{output_dir}/$pref\_by_clone.txt | grep ';$config->{lineage2};' | wc -l");
    my $valid_int_count = $self->_run_cmd("wc -l $self->{output_dir}/$pref\_by_trace.txt");
    chomp $valid_count;
    chomp $valid_int_count;
    return {valid_clones => $valid_count,
            valid_traces => $valid_int_count
    };
}

=head2 splitBam

 Title   : splitBam
 Usage   : $lgtseek->splitBam(({'input' => $file})
 Function: Split a bam file into smaller chunks.
 Returns : A list of the bam files split up
 Args    : An object with the input bam files.

=cut
sub splitBam {
    my ($self,$config) = @_;

    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    my $seqs_per_file = $config->{seqs_per_file} ? $config->{seqs_per_file} : 10000;

    my $samtools = $self->{samtools_bin};

    # Parse out the pieces of the filename
    my($fn,$path,$suffix) = fileparse($config->{input},'.bam');

    # Pull out the header
    my $header = $self->_run_cmd("$samtools view -H $config->{input}");

    # Open the input file
    open(IN, "-|", "$samtools view $config->{input}");

    # Open the first output file
    my $count = 0;
    my $ofile = "$output_dir/$fn\_$count.bam";
    open(my $ofh, "| samtools view -S -b -o $ofile -");
    my @outfiles = ($ofile);

    # Print the header
    print $ofh $header;
    
    my $i = 0;
    while(my $line = <IN>) {
        my @fields = split(/\t/,$line);
        my $flag = $self->_parseFlag($fields[1]);

        # Strip out the XA tag
        $line =~ s/\s+XA:Z:\S+//;

        # Make sure we don't accidentally split read pairs.
        if($flag->{'last'} && $i >= $seqs_per_file) {

            # Print out the line
            print $ofh $line;

            # Close the old file
            close $ofh;
            
            # Open the new file
            $count ++;
            $ofile = "$output_dir/$fn\_$count.bam";
            open($ofh, "| samtools view -S -b -o $ofile -");
            push(@outfiles,$ofile);

            # Print the header to the new file
            print $ofh $header;
            
            # Reset the counter
            $i = 0;
        }

        # If we haven't filled the current file yet just keep printing.
        else {
            print $ofh $line;
        }

        # Increment the counter
        $i++;
    }
    close $ofh;
    my $files = \@outfiles;
    return $files; 
}

sub _dec2bin {
    my $self = shift;
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub _parseFlag {
    my $self = shift;
    my $int = shift;
    my $rawbin = $self->_dec2bin($int);
    my $rev = scalar $rawbin;
    if($rev eq $rawbin) {
    #    print "ERROR $rev $rawbin\n";
    }
    my $bin = sprintf("%011d", $rev);
    my $final_bin = reverse $bin;
    return {
        'paired' => substr($final_bin, 0, 1),
        'proper' => substr($final_bin, 1, 1),
        'qunmapped' => substr($final_bin, 2, 1),
        'munmapped' => substr($final_bin, 3, 1),
        'qrev' => substr($final_bin, 4, 1),
        'mrev' => substr($final_bin, 5, 1),
        'first' => substr($final_bin, 6, 1),
        'last' => substr($final_bin, 7, 1),
        'secondary' => substr($final_bin, 8, 1),
        'failqual' => substr($final_bin, 9, 1),
        'pcrdup' => substr($final_bin, 10, 1)
    };
}
1;
