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

# Dependencies
use GiTaxon;
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

    if(!$self->{ergatis_bin} || $self->{prinseq_bin}) {
        die "Must provide an ergatis_bin and prinseq_bin parameter to run prinseq filtering\n";
    }

    my $retval;
    if($self->{paired_end}) {
        $retval = $self->_prinseqFilterPaired($config->{bam_file});
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
    my($self,$bam_file) = @_;

    my($name,$path,$suff) = fileparse($bam_file,'.bam');

    my $bin = $self->{ergatis_bin};
    my $prinseq_bin = $self->{prinseq_bin};
    my $samtools = $self->{samtools_bin};

    # Generate concatenated fastq files for prinseq derep filtering
    my $cmd = "perl $bin/sam2fasta.pl --input=$bam_file --fastq=1 --combine_mates=1 --output_file=$path/$name\_combined.fastq";
    $self->_run_cmd($cmd);
    # Run prinseq for dereplication
    my $cmd = "perl $prinseq_bin --fastq=$path/$name\_combined.fastq --out_good=$path/$name\_derep_good --out_bad=$path/$name\_derep_bad -derep 14";
    $self->_run_cmd($cmd);

    # Pull out bad ids
    my $cmd = "perl -e 'while(<>){s/\@//;print;<>;<>;<>;}' $path/$name\_derep_bad.fastq > $path/$name\_derep_bad_ids.out";
    print STDERR "$cmd\n";
    $self->_run_cmd($cmd);

    # Generate single-read fastq for low complexity filtering
    my $cmd = "perl $bin/sam2fasta.pl --input=$bam_file --fastq=1 --combine_mates=0 --paired=1 --output_file=$path/$name.fastq";
    $self->_run_cmd($cmd);

    # Run prinseq for low complexity filtering
    my $cmd = "perl /mnt/prinseq-lite2.pl --fastq=$path/$name\_1.fastq --out_good=$path/$name\_lc_1_good --out_bad=$path/$name\_lc_1_bad -lc_method dust -lc_threshold 7";
    $self->_run_cmd($cmd);

    # Pull out bad ids
    my $cmd = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $path/$name\_lc_1_bad.fastq > $path/$name\_lc_1_bad_ids.out";
    $self->_run_cmd($cmd);


    # Run prinseq for low complexity filtering
    my $cmd = "perl /mnt/prinseq-lite2.pl --fastq=$path/$name\_2.fastq --out_good=$path/$name\_lc_2_good --out_bad=$path/$name\_lc_2_bad -lc_method dust -lc_threshold 7";
    $self->_run_cmd($cmd);

    # Pull out bad ids
    my $cmd = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $path/$name\_lc_2_bad.fastq > $path/$name\_lc_2_bad_ids.out";
    $self->_run_cmd($cmd);

    # Merge bad ids from derep and lc filtering
    my $cmd = "cat $path/$name\_derep_bad_ids.out $path/$name\_lc_1_bad_ids.out $path/$name\_lc_2_bad_ids.out | sort -u > $path/$name\_bad_ids.out";
   $self->_run_cmd($cmd);

    # Filter sam file to remove bad ids

    my $cmd = "perl $bin/filter_sam_from_prinseq.pl --sam_file=$bam_file --bad_list=$path/$name\_bad_ids.out --out_file=$path/$name\_filtered.sam";
    $self->_run_cmd($cmd);

    my $cmd = "$samtools view -S -b $path/$name\_filtered.sam > $path/$name\_filtered.bam";
    $self->_run_cmd($cmd);

    $self->_run_cmd($cmd);
#    `rm $path/$name\_filtered.sam`; 
    return "$path/$name\_filtered.bam";

}


=head2 run_cmd

 Title   : _run_cmd
 Usage   : *PRIVATE*
 Function: Run a unix command and fail if something goes wrong
 Returns : void
 Args    : Command to run

=cut
sub _run_cmd {

    my $cmd = shift;

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

    my $cmd_string = "$self->{aspera_path}/bin/ascp -QTd -l$self->{aspera_rate} -i $self->{aspera_path}/etc/asperaweb_id_dsa.putty $self->{aspera_user}\@$self->{aspera_host}:$path_to_file $self->{output_dir} -L $output_dir -o Overwrite=diff 2>&1";

    # Doing this echo y to ensure we accept any certs.
    my $out = `echo y | $cmd_string`;

    # We can actually exit non-0 and still succeed if the 
    if($out =~ /Error/)  {
        print STDERR "Had a problem downloading $self->{aspera_host}:$path_to_file to $output_dir\n";
        print STDERR "$cmd_string";
        exit(1);
    }
    else {
        print STDERR "$? $out Successfully downloaded $self->{aspera_host}:$path_to_file to $output_dir\n";
    }
    my @files = `find $output_dir -name '*.sra'`;
    return \@files;
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


    # Need to pull the version of the sratoolkit to determine if we need the --split-3 parameter.
    my $ret = `$self->{sratoolkit_path}/fastq-dump -V`;
    $ret =~ /fastq-dump : ([\d.]+)/;
    my $version = version->parse($1);
    my $cutoff_version = version->parse('2.1.0');

    my $fastqdump_bin = "$self->{sratoolkit_path}/fastq-dump";

    if($version > $cutoff_version) {

        $fastqdump_bin .= " --split-3 ";
    }
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if(!$self->{output_dir}) {
        die "Need to specify an output_dir in order to download from the sra\n";
    }

    my ($name,$path,$suff) = fileparse($config->{sra_file});
    my $cmd = "$fastqdump_bin -O $self->{output_dir} $config->{sra_file}";

    $self->_run_cmd($cmd);

    my $res = $self->_run_cmd("find $self->{output_dir} -name '$name*.fastq'");
    my @files = split(/\n/,$res);

    my $retval = {
        'files' => \@files,
        'path' => $self->{output_dir},
        'basename' => $name,
        'paired_end' => 0
    };
    if($files[0] =~ /_1.fastq/) {
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

    $self->{bin_dir} = $config->{bin_dir} ? $config->{bin_dir} : $self->{bin_dir};
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    # Check for a bwa path. If we don't have one we'll just hope it's in our global path.
    $self->{bwa_path} = $config->{bwa_path} ? $config->{bwa_path} : $self->{bwa_path};
    $self->{bwa_path} = $self->{bwa_path} ? $self->{bwa_path} : 'bwa';


    # Build the command string;
    my @cmd = ("$self->{bin_dir}/lgt_bwa --bwa_path=$self->{bwa_path} --bin_dir=$self->{bin_dir}");

    my $suff = '.sam';
    if($config->{output_bam}) {
        $suff = '.bam';
        push(@cmd, '--output_bam=1');
    }

    my $basename = $config->{input_base};

    # Handle making up the lgt_bwa command with a bam file
    if($config->{input_bam}) {
        my ($name,$path,$suff) = fileparse($config->{input_bam},'.bam');
        $basename = $name;
        
        push(@cmd,"--input_bam=$config->{input_bam}");
    }
    elsif($config->{input_dir} && $config->{input_base}) {
        push(@cmd,"--input_dir=$config->{input_base} --input_base=$config->{input_base}");
    }
    else {
        die "Must provide either a value to either input_bam or to input_base and input_dir\n";
    }

    if($config->{reference}) {
        push(@cmd,"--ref_file=$config->{reference}");
    }
    elsif($config->{reference_list}) {
        push(@cmd,"--ref_file_list=$config->{reference_list}");  
    }
    else {
        die "Must provide a value for either reference or reference_list to run bwa\n";
    }

    push(@cmd, $self->{output_dir});
    push(@cmd, $config->{other_opts});

    my $cmd_string = join(' ',@cmd);

    $self->_run_cmd($cmd_string);
    
    my @files = $self->_run_cmd("find $self->{output_dir} -name '*$basename$suff'");
    
    return \@files;
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

    # Do we have both donor and host bams?
    if($config->{donor_bams} && $config->{host_bams}) {
        $self->_bwaPostProcessDonorHostPaired($config);
    }
}


sub _bwaPostProcessDonorHostPaired {
    my ($self, $config) = @_;

    my $samtools = $self->{samtools_bin};

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
        'MM_UU' => 'microbiome'
    };

    my $prefix = $config->{output_prefix} ? "$config->{output_prefix}_" : '';
    # Here are a bunch of file handles we'll use later.
    print STDERR "$config->{output_dir}/".$prefix."lgt_donor.bam\n";
    open(my $lgtd,  "| $samtools view -S -b -o $config->{output_dir}/".$prefix."lgt_donor.bam -") or die "Unable to open\n";
    open(my $lgth, "| $samtools view -S -b -o $config->{output_dir}/".$prefix."lgt_host.bam -") or die "Unable to open\n";
    open(my $int_site_donor_d, "| $samtools view -S -b -o $config->{output_dir}/".$prefix."integration_site_donor_donor.bam -") or die "Unable to open\n";
    open(my $int_site_donor_h, "| $samtools view -S -b -o $config->{output_dir}/".$prefix."integration_site_donor_host.bam -") or die "Unable to open\n";
    open(my $microbiome_donor,"| $samtools view -S -b -o $config->{output_dir}/".$prefix."microbiome.bam -") or die "Unable to open\n";

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
        push(@donor_head, `$samtools view -H $_`);
        open(my $fh, "-|", "$samtools view $_");
        push(@donor_fh,$fh);
    } @{$config->{donor_bams}};
    
    
    
    # Open all the host files
    map {
        print STDERR "Opening $_\n";
        push(@host_head, `$samtools view -H $_`);
        open(my $fh, "-|", "$samtools view $_");
        push(@host_fh,$fh);
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

    my $class_counts = {
        'lgt' => 0,
        'integration_site_host' => 0,
        'integration_site_donor' => 0,
        'microbiome' => 0
    };

    while($more_lines) {
        
        my @donor_lines;
        my @host_lines;
        
        my $hr1_class;
        my $hr2_class;
        my $hr1_line;
        my $hr2_line;

        # First establish the class of the host reads.
        foreach my $hfh (@host_fh) {
            my $hr1 = <$hfh>;
            my $hr2 = <$hfh>;            
            
            # Should check if these ended at the same time?
            if(!$hr1 || !$hr2 || $more_lines == 0) {
                $more_lines = 0;
                last;
            }
            my $hr1_flag = $self->_parseFlag((split(/\t/,$hr1))[1]);
 #           my $hr2_flag = $self->_parseFlag((split(/\t/,$hr2))[1]);
            if(!$hr1_flag->{'qunmapped'}) {
                $hr1_line = $hr1;
                $hr1_class = 'M';
            }
            elsif(!$hr1_class) {
                $hr1_line = $hr1;
                $hr1_class = 'U';
            }
            if(!$hr1_flag->{'munmapped'}) {
                $hr2_line = $hr2;
                $hr2_class = 'M';
            }
            elsif(!$hr2_class) {
                $hr2_line = $hr2;
                $hr2_class = 'U';
            }
        }

        my $dr1_class;
        my $dr2_class;
        my $dr1_line;
        my $dr2_line;

        # Next establish the class of the donor read
        foreach my $dfh (@donor_fh) {
            my $dr1 = <$dfh>;
            my $dr2 = <$dfh>;
#                print STDERR "Processing $hr1$hr2$dr1$dr2";
            # Should check if these ended at the same time?
            if(!$dr1 || !$dr2) {
                $more_lines = 0;
                last;
            }

            my $dr1_flag = $self->_parseFlag((split(/\t/,$dr1))[1]);
#            my $dr2_flag = $self->_parseFlag((split(/\t/,$dr2))[1]);
            if(!$dr1_flag->{'qunmapped'}) {
                $dr1_line = $dr1;
                $dr1_class = 'M';
            }
            elsif(!$dr1_class) {
                $dr1_line = $dr1;
                $dr1_class = 'U';
            }
            if(!$dr1_flag->{'munmapped'}) {
                $dr2_line = $dr2;
                $dr2_class = 'M';
            }            
            elsif(!$dr2_class) {
                $dr2_line = $dr2;
                $dr2_class = 'U';
            }
        }

        my $dclass = "$dr1_class$dr2_class";
        my $hclass = "$hr1_class$hr2_class";
        my $paired_class = "$dclass\_$hclass";
        
        if($class_to_file->{$classes_both->{$paired_class}."_donor"}) {
            print STDERR "Printing to ".$classes_both->{$paired_class}."_donor $paired_class\n$dr1_line$dr2_line";
            print {$class_to_file->{$classes_both->{$paired_class}."_donor"}} "$dr1_line$dr2_line";
        }
        if($class_to_file->{$classes_both->{$paired_class}."_host"}) {
            print {$class_to_file->{$classes_both->{$paired_class}."_host"}} "$hr1_line$hr2_line";
        }

        if($classes_both->{$paired_class} eq 'lgt') {
#            print STDERR "Processing $hr1_line$hr2_line$dr1_line$dr2_line";
        }
        if($classes_both->{$paired_class}) {
            $class_counts->{$classes_both->{$paired_class}}++;
            my @lines;
#            map {push(@lines, "$_: $class_counts->{$_}")} keys %$class_counts;
#            map {print  "\r$_: $class_counts->{$_}"} keys %$class_counts;
#                   print "\r".join(' ',@lines);
                
                #print STDERR "Line $line_num donor class: $dclass\nHost class: $hclass\nCombined class: $paired_class\n";
        }
        $line_num ++;
    }

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
