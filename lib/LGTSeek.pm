
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

=head2 LGTSeek.pm

 Title   : LGTSeek.pm
 Usage   : Suite of subroutines used to identify LGT
 Routines:
        new                 : Create lgtseek object
        downloadSRA         : Download SRA file
        dumpFastq           : Convert sra file 2 fastq
        downloadCGHub       : Download bam from GC-Hub
        decrypt             : 
        prelim_filter       : Filter for potential LGT reads from human mapped bam
        runBWA              : Execute BWA mapping
        bwaPostProcess      : 
        prinseqFilterBam    :
        filter_bam_by_ids   :
        sam2fasta           :
        splitBam            :
        blast2lca           :
        bestBlast2          :
        runLgtFinder        :
        getGiTaxon          : 
        mpileup             :
        empty_chk           :
        fail                :
        _parseFlag          :
        _run_cmd            :
=cut

package LGTSeek;
## use warnings;
use strict;
use version;
use Carp;
$Carp::MaxArgLen = 0;
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
    my ( $class, $args ) = @_;

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
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &getGiTaxon: Start ========\n"; }

    # If we already have a gitaxon object we'll just return it.
    if ( !$self->{gitaxon} ) {

        # Apply any config options that came over.
        if ($config) {
            map { $self->{$_} = $config->{$_} } keys %$config;
        }

        # Create the object.
        $self->{gitaxon} = GiTaxon->new(
            {   'taxon_dir'  => $self->{taxon_dir},
                'chunk_size' => 10000,
                'idx_dir'    => $self->{taxon_idx_dir},
                'host'       => $self->{taxon_host},
                'type'       => 'nucleotide'
            }
        );
    }
    if ( $self->{verbose} ) { print STDERR "======== &getGiTaxon: Finished ========\n"; }
    return $self->{gitaxon};

}

=head2 &prinseqFilterBam

 Title   : prinseqFilterBam
 Usage   : my $filteredBam = $LGTSeek->prinseqFilterBam({'input_bam' => '/path/to/file.bam'...})
 Function: Prinseq filter a bam file
 Returns : A Hash $ref->{bam} = path to the filtered bam. $ref->{count} = # of reads passing filtering.
 Args    : 
        input_bam    => /path/to/file.bam
        output_dir   => /path/for/output.bam
        overwrite    => <0|1> [0] 1= Overwrite output if it is found already. 
        prinseq_bin  =>
        samtools_bin =>
        ergatis_bin  =>
        bin_dir      =>
=cut

sub prinseqFilterBam {
    my ( $self, $config ) = @_;
    if ( !$config->{input_bam} ) { die "Must pass &prinseqFilterBam and input =>\n" }
    if ( $self->empty_chk( { input => $config->{input_bam} } ) == 1 ) {
        print STDERR "Warning: &prinseqFilterBam input: $config->{input} is empty.\n";
    }
    if ( $self->{verbose} ) { print STDERR "======== &prinseqFilterBam: Start ========\n"; }

    # Override if it is provided
    $self->{prinseq_bin}  = $config->{prinseq_bin} ? $config->{prinseq_bin} : $self->{prinseq_bin};
    $self->{samtools_bin} = $self->{samtools_bin}  ? $self->{samtools_bin}  : 'samtools';
    my $overwrite = $config->{overwrite} ? $config->{overwrite} : 0;
    if ( $config->{output_dir} ) {
        $self->_run_cmd("mkdir -p $config->{output_dir}");
    }

    if ( !$self->{ergatis_bin} || !$self->{prinseq_bin} ) {
        $self->fail(
            "Must provide an ergatis_bin and prinseq_bin parameter to run prinseq filtering $self->{prinseq_bin} $self->{ergatis_bin}\n"
        );
    }

    my $retval;
    if ( $self->{paired_end} ) {
        $retval = $self->_prinseqFilterPaired( $config->{input_bam}, $config->{output_dir}, $overwrite );
    }
    else {
        $self->fail("Single end is currently not implemented\n");
    }
    if ( $self->{verbose} ) { print STDERR "======== &prinseqFilterBam: Finished ========\n"; }
    return $retval;
}

=head2 &_prinseqFilterPaired

 Title   : _prinseqFilterPaired
 Usage   : *PRIVATE*
 Function: Prinseq filter a bam paired end file
 Returns : Path to the filtered bam file
 Args    : @_=($self,<bam_to_filter>,<output_dir>)

=cut

sub _prinseqFilterPaired {
    my ( $self, $bam_file, $output_dir, $overwrite ) = @_;

    my ( $name, $path, $suff ) = fileparse( $bam_file, ".bam" );

    $output_dir = $output_dir ? $output_dir : $path;

    my $tmp_dir = "$output_dir/tmp/";

    $self->_run_cmd("mkdir -p $output_dir");
    $self->_run_cmd("mkdir -p $tmp_dir");

    my $bin         = $self->{bin_dir};
    my $prinseq_bin = $self->{prinseq_bin};

    my $samtools = $self->{samtools_bin};

    if ( -e "$output_dir/$name\_bad_ids.out" && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR "Already found the output for &prinseqFilter: $output_dir/$name\_bad_ids.out";
        }
        my $filtered = $self->filter_bam_by_ids(
            {   input_bam => $bam_file,
                bad_list  => "$output_dir/$name\_bad_ids.out",
            }
        );
        return $filtered;
    }
    else {
        # Generate concatenated fastq files for prinseq derep filtering
        ## Need to incorporate sam2fasta.pm KBS 01.07.14
        if ( $self->{verbose} ) { print STDERR "======== Deduplication Filtering========\n"; }
        my $cmd
            = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$bam_file --fastq=1 --combine_mates=1 --output_file=$tmp_dir/$name\_combined.fastq";

        # print STDERR "$cmd\n";
        $self->_run_cmd($cmd);

        # Run prinseq for dereplication
        $cmd
            = "perl $prinseq_bin --fastq=$tmp_dir/$name\_combined.fastq --out_good=$tmp_dir/$name\_derep_good --out_bad=$tmp_dir/$name\_derep_bad -derep 14";

        # print STDERR "$cmd\n";
        $self->_run_cmd($cmd);

        # Pull out bad ids
        if ( $self->{verbose} ) { print STDERR "======== Pull DeDup Bad ID's ========\n"; }
        $cmd
            = "perl -e 'while(<>){s/\@//;print;<>;<>;<>;}' $tmp_dir/$name\_derep_bad.fastq > $tmp_dir/$name\_derep_bad_ids.out";

        # print STDERR "$cmd\n";
        $self->_run_cmd($cmd);

        # Generate single-read fastq for low complexity filtering
        if ( $self->{verbose} ) { print STDERR "======== Low Complexity Filter ========\n"; }
        $cmd
            = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$bam_file --fastq=1 --combine_mates=0 --paired=1 --output_file=$tmp_dir/$name.fastq";
        $self->_run_cmd($cmd);

        # Run prinseq for low complexity filtering
        $cmd
            = "perl $prinseq_bin --fastq=$tmp_dir/$name\_1.fastq --out_good=$tmp_dir/$name\_lc_1_good --out_bad=$tmp_dir/$name\_lc_1_bad -lc_method dust -lc_threshold 7";
        $self->_run_cmd($cmd);

        if ( -e "$tmp_dir/$name\_lc_1_bad.fastq" ) {

            # Pull out bad ids
            if ( $self->{verbose} ) { print STDERR "======== Pull Low-Cmplx-1 Bad ID's ========\n"; }
            my $cmd
                = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_1_bad.fastq > $tmp_dir/$name\_lc_1_bad_ids.out";
            $self->_run_cmd($cmd);
        }
        else {
            if ( $self->{verbose} ) { print STDERR "Didn't find any low complexity sequences in read 1\n"; }
            $self->_run_cmd("touch $tmp_dir/$name\_lc_1_bad_ids.out");
        }

        # Run prinseq for low complexity filtering
        $cmd
            = "perl $prinseq_bin --fastq=$tmp_dir/$name\_2.fastq --out_good=$tmp_dir/$name\_lc_2_good --out_bad=$tmp_dir/$name\_lc_2_bad -lc_method dust -lc_threshold 7";
        $self->_run_cmd($cmd);

        # Pull out bad ids
        if ( -e "$tmp_dir/$name\_lc_2_bad.fastq" ) {
            if ( $self->{verbose} ) { print STDERR "======== Pull Low-Cmplx-2 Bad ID's ========\n"; }
            my $cmd
                = "perl -e 'while(<>){s/\@//;s/\_\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_2_bad.fastq > $tmp_dir/$name\_lc_2_bad_ids.out";
            $self->_run_cmd($cmd);
        }
        else {
            if ( $self->{verbose} ) { print STDERR "Didn't find any low complexity sequences in read 2\n"; }
            $self->_run_cmd("touch $tmp_dir/$name\_lc_2_bad_ids.out");
        }

        # Merge bad ids from derep and lc filtering
        $cmd
            = "cat $tmp_dir/$name\_derep_bad_ids.out $tmp_dir/$name\_lc_1_bad_ids.out $tmp_dir/$name\_lc_2_bad_ids.out | sort -u > $output_dir/$name\_prinseq-bad-ids.out";

        $self->_run_cmd($cmd);

        # Filter bam file to remove bad ids

        my $filtered = $self->filter_bam_by_ids(
            {   input_bam  => $bam_file,
                output_dir => $output_dir,
                bad_list   => "$output_dir/$name\_prinseq-bad-ids.out",
            }
        );

        # Cleanup /tmp/ folder w/ intermediate files
        $self->_run_cmd("rm -rf $tmp_dir");
        return $filtered;
    }
    ## Incorporated filter_sam_from_pringseq.pl into &filter_bam_by_ids
#my $cmd = "perl $bin/filter_sam_from_prinseq.pl --sam_file=$bam_file --bad_list=$output_dir/$name\_bad_ids.out --out_file=$output_dir/$name\_filtered.sam";
#$self->_run_cmd($cmd);

    #my $cmd = "$samtools view -S -b $output_dir/$name\_filtered.sam > $output_dir/$name\_filtered.bam";
    #$self->_run_cmd($cmd);

    #my $count = $self->_run_cmd("$samtools view $output_dir/$name\_filtered.bam | cut -f1 | uniq | wc -l");
    #chomp $count;

    # Blitz the sam file
    #$self->_run_cmd("rm $output_dir/$name\_filtered.sam");

#    my $cmd = "perl $bin/sam2fasta.pl --input=$output_dir/$name\_filtered.bam --combine_mates=0 --paired=1 --output_file=$output_dir/$name.fastq";
#    $self->_run_cmd($cmd);
#    `rm $output_dir/$name\_filtered.sam`;
    ## Return ->{count} and ->{file}

}

=head2 &sam2Fasta
## This needs to be updated to incorporate the sam2fasta.pm into LGTSeek.pm OR use sam2fastaconverter.pm
 Title   : sam2Fasta
 Usage   : my $fastas = $LGTSeek->sam2Fasta({'input' => '/path/to/file.bam'...})
 Function: Convert a bam/sam file to a fasta file
 Returns : a list of fasta/fastq files
 Args    : input => sam or bam file to convert to fasta
           output_dir => directory for output

=cut

sub sam2Fasta {
    my ( $self, $config ) = @_;
    my $bin = $self->{ergatis_bin};

    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) { print STDERR "======== &sam2Fasta: Start ========\n"; }

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");

    my $outfile;
    my ( $name, $path, $suff ) = fileparse( $config->{input}, (qr/_filtered.bam$||.bam$||.sam$||.sam.gz$/) );
    my $cmd = "perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$config->{input}";
    if ( $config->{fastq} ) {
        $outfile = "$output_dir/$name.fastq";
        $cmd .= " --fastq=1 --output_file=$outfile";
    }
    else {
        $outfile = "$output_dir/$name.fasta";
        $cmd .= " --fastq=0 --output_file=$outfile";
    }
    if ( $config->{combine_mates} ) {
        $cmd .= " --combine_mates=0";
    }
    if ( $config->{paired} || $self->{paired_end} ) {
        $cmd .= " --paired=1";
    }

    $self->_run_cmd($cmd);
    if ( $self->{verbose} ) { print STDERR "======== &sam2Fasta: Finished ========\n"; }
    return $outfile;
}

=head2 &_run_cmd

 Title   : _run_cmd
 Usage   : *PRIVATE*
 Function: Run a unix command and fail if something goes wrong
 Returns : void
 Args    : Command to run

=cut

sub _run_cmd {

    my ( $self, $cmd ) = @_;

    if ( $self->{verbose} ) { print STDERR "CMD: $cmd\n"; }
    my $res = `$cmd`;
    if ($?) {
        print STDERR "$cmd\n\n$?";
        $self->fail("$cmd died with message:\n$res\n\n");
    }
    return $res;
}

=head2 &downloadSRA

 Title   : downloadSRA
 Usage   : $lgtseek->downloadSRA(({'experiment_id'} = 'SRX01234'})
 Function: Download sra files from the sequence read archive
 Returns : A list of the downloaded file paths
 Args    : 

=cut

sub downloadSRA {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &downloadSRA: Start ========\n"; }

    # Check for all the aspera related options
    $self->{aspera_host}   = $config->{aspera_host}   ? $config->{aspera_host}   : $self->{aspera_host};
    $self->{aspera_user}   = $config->{aspera_user}   ? $config->{aspera_user}   : $self->{aspera_user};
    $self->{aspera_params} = $config->{aspera_params} ? $config->{aspera_params} : $self->{aspera_params};
    $self->{aspera_path}   = $config->{aspera_path}   ? $config->{aspera_path}   : $self->{aspera_path};
    $self->{aspera_rate}   = $config->{aspera_rate}   ? $config->{aspera_rate}   : $self->{aspera_rate};
    $self->{output_dir}    = $config->{output_dir}    ? $config->{output_dir}    : $self->{output_dir};
    $self->{aspera_path}   = $self->{aspera_path}     ? $self->{aspera_path}     : '~/.aspera/connect/';

    my $retry_attempts = $config->{retry_attempts} ? $config->{retry_attempts} : 10;
    if ( !$self->{aspera_host} ) {
        $self->{aspera_host} = 'ftp-private.ncbi.nlm.nih.gov';
    }
    if ( !$self->{aspera_user} ) {
        $self->{aspera_user} = 'anonftp';
    }
    if ( !$self->{aspera_rate} ) {
        $self->{aspera_rate} = '200M';
    }

    if ( !$self->{aspera_path} ) {
        $self->fail("Need to specify an aspera_path (where is aspera installed) in order to download from the sra\n");
    }
    if ( !$self->{output_dir} ) {
        $self->fail("Need to specify an output_dir in order to download from the sra\n");
    }

    my $prefix;
    my $thousand;
    my $output_dir = $self->{output_dir};
    my $path_to_file;

    # We can pass an experiment_id, run_id or a full path to download
    if ( $config->{experiment_id} ) {
        my $exp_id = $config->{experiment_id};
        $exp_id =~ /^((.{3}).{3}).*/;
        $prefix       = $2;
        $thousand     = $1;
        $path_to_file = "/sra/sra-instant/reads/ByExp/litesra/$prefix/$thousand/$exp_id";
        $output_dir   = "$output_dir/$prefix/$thousand/";
    }
    if ( $config->{run_id} ) {
        $config->{run_id} =~ /^((.{3}).{3}).*/;
        $prefix       = $2;
        $thousand     = $1;
        $path_to_file = "/sra/sra-instant/reads/ByRun/litesra/$prefix/$thousand/$config->{run_id}";
        $output_dir   = "$output_dir/$prefix/$thousand/";
    }
    elsif ( $config->{path} ) {
        $path_to_file = $config->{path};
    }

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");

    my $cmd_string
        = "$self->{aspera_path}/bin/ascp -QTd -l$self->{aspera_rate} -i $self->{aspera_path}/etc/asperaweb_id_dsa.putty $self->{aspera_user}\@$self->{aspera_host}:$path_to_file $output_dir -L $output_dir -o Overwrite=diff 2>&1";

    #Retry the download several times just incase.
    my $retry = 1;

    my $retries = 0;
    while ($retry) {

        # Doing this echo y to ensure we accept any certs.
        my $out = $self->_run_cmd("echo y | $cmd_string");

        # We can actually exit non-0 and still succeed if the
        if ( $out =~ /Error/ ) {
            print STDERR "Had a problem downloading $self->{aspera_host}:$path_to_file to $output_dir\n";
            print STDERR "$cmd_string";
            if ( $retries < $retry_attempts ) {
                $retries++;
                sleep $retries * 2;    # Sleep for 2 seconds per retry.
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
    if ( $self->{verbose} ) { print STDERR "======== &downloadSRA: Finished ========\n"; }
    return \@files;
}

=head2 &dumpFastq

 Title   : dumpFastq
 Usage   : $lgtseek->dumpFastq(({'sra_file' => 'SRX01234'})
 Function: Run the sratoolkit program dump-fastq
 Returns : An object with a path and basename of the output as well as a list of output files
 Args    : An object with element 'sra_file' and optionally the path to the sratoolkit install

=cut

sub dumpFastq {
    my ( $self, $config ) = @_;

    $self->{sratoolkit_path} = $config->{sratoolkit_path} ? $config->{sratoolkit_path} : $self->{sratoolkit_path};
    if ( $self->{verbose} ) { print STDERR "======== &dumpFastq: Start ========\n"; }

    # If we don't have a path provided we'll hope it's in our path.
    my $fastqdump_bin = $self->{sratoolkit_path} ? "$self->{sratoolkit_path}/fastq-dump" : "fastq-dump";

    $config->{sra_file} =~ s/\/\//\//g;

    # Need to pull the version of the sratoolkit to determine if we need the --split-3 parameter.
    my $ret = `$fastqdump_bin -V`;
    my $version;
    my $cutoff_version;
    if ( $ret =~ /fastq-dump : ([\d.]+)/ ) {
        $version        = version->parse($1);
        $cutoff_version = version->parse('2.1.0');
    }
    else {
        $self->fail("$? $ret $fastqdump_bin\n");
    }
    if ( $version > $cutoff_version ) {

        $fastqdump_bin .= " --split-3 ";
    }
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( !$self->{output_dir} ) {
        $self->fail("Need to specify an output_dir in order to download from the sra\n");
    }

    my ( $name, $path, $suff ) = fileparse( $config->{sra_file}, ".sra" );
    chomp $name;
    my $res = $self->_run_cmd("find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
    my @files = split( /\n/, $res );

    if ( !@files && !$config->{overwrite} ) {
        my $cmd = "$fastqdump_bin -O $self->{output_dir} $config->{sra_file}";
        $self->_run_cmd($cmd);
        my $res = $self->_run_cmd("find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
        @files = split( /\n/, $res );
    }

    if ( $self->{verbose} ) { print STDERR "@files\n"; }

    my $retval = {
        'files'      => \@files,
        'path'       => $self->{output_dir},
        'basename'   => $name,
        'paired_end' => 0
    };
    if ( $files[0] =~ /_\d.fastq/ ) {
        $retval->{paired_end} = 1;
    }
    if ( $self->{verbose} ) { print STDERR "======== &dumpFastq: Finished ========\n"; }
    return $retval;
}

=head2 downloadCGHub

 Title   : downloadCGHub
 Usage   : $lgtseek->downloadCGHub(({'analysis_id'} = '00007994-abeb-4b16-a6ad-7230300a29e9'})
 Function: Download TCGA bam files from the CGHub
 Returns : A list of the downloaded file paths
 Args    : 

=cut

sub downloadCGHub {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &downloadCGHub: Start ========\n"; }

    # Check for all the genetorrent related options
    $self->{output_dir}       = $config->{output_dir}     ? $config->{output_dir}     : $self->{output_dir};
    $self->{genetorrent_path} = $self->{genetorrent_path} ? $self->{genetorrent_path} : '/opt/genetorrent/bin';
    $self->{cghub_key}        = $config->{cghub_key}      ? $config->{cghub_key}      : $self->{cghub_key};

    my $retry_attempts = $config->{retry_attempts} ? $config->{retry_attempts} : 10;
    if ( !$self->{cghub_key} ) {

        #        $self->{aspera_user} = 'anonftp';
        die "Need to specify the path to the cghub_key\n";
    }

    if ( !$self->{genetorrent_path} ) {
        die
            "Need to specify an genetorrent_path (where is genetorrent/gtdownload installed) in order to download from CGHub\n";
    }
    if ( !$self->{output_dir} ) {
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
    if ( $self->{verbose} ) { print STDERR "======== &downloadCGHub: Finished ========\n"; }
    return $files;
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
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &runBWA: Start ========\n"; }
    $self->{ergatis_bin} = $config->{ergatis_bin} ? $config->{ergatis_bin} : $self->{ergatis_bin};
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    # Check for a bwa path. If we don't have one we'll just hope it's in our global path.
    $self->{bwa_path} = $config->{bwa_path} ? $config->{bwa_path} : $self->{bwa_path};
    $self->{bwa_path} = $self->{bwa_path}   ? $self->{bwa_path}   : 'bwa';
    my $num_threads = defined $config->{threads} ? $config->{threads} : $self->{threads};
    $self->_run_cmd("mkdir -p $output_dir");

    my $conf = {
        num_aligns => 3,
        bwa_path   => $self->{bwa_path},
        output_dir => $output_dir,
        threads    => $num_threads,
    };

    # Build the command string;
    my @cmd = ("$self->{ergatis_dir}/lgt_bwa --num_aligns=3 --bwa_path=$self->{bwa_path}");

    my $suff = '.sam';
    if ( $config->{output_bam} ) {
        $suff = '.bam';
        $conf->{output_bam} = 1;
    }

    my $basename = $config->{input_base};

    # Handle making up the lgt_bwa command with a bam file
    if ( $config->{input_bam} ) {
        if ( $self->empty_chk( { input => $config->{input_bam} } ) == 1 ) {
            print STDERR "Error: &runBWA input: $config->{input_bam} is empty.\n";
            return $config->{input_bam};
        }
        my ( $name, $path, $suff ) = fileparse( $config->{input_bam}, ( "_prelim.bam", ".bam" ) );
        $basename           = $name;
        $conf->{input_base} = $basename;
        $conf->{input_bam}  = $config->{input_bam};
    }
    elsif ( $config->{input_dir} && $config->{input_base} ) {
        $conf->{input_dir}  = $config->{input_dir};
        $conf->{input_base} = $config->{input_base};
    }
    else {
        die "Must provide either a value to either input_bam or to input_base and input_dir\n";
    }

    my $pre = '';
    if ( $config->{reference} ) {
        my ( $name, $dir, $suff ) = fileparse( $config->{reference}, qr/\.[^\.]+/ );
        $pre = "$name\_";
        $conf->{ref_file} = $config->{reference};
    }
    elsif ( $config->{reference_list} ) {
        $conf->{ref_file_list} = $config->{reference_list};
    }
    else {
        die "Must provide a value for either reference or reference_list to run bwa\n";
    }

    $conf->{overwrite} = $config->{overwrite};
    map { $conf->{$_} = $config->{other_opts}->{$_}; } keys %{ $config->{other_opts} };
    $conf->{run_lca}     = $config->{run_lca};
    $conf->{lgtseek}     = $self;
    $conf->{cleanup_sai} = $config->{cleanup_sai};
    $conf->{out_file}    = $config->{out_file};
    LGTbwa::runBWA($conf);

    # Maybe should check if this is valid.
    if ( $config->{run_lca} ) {

    }
    else {
        my @files = split( /\n/, $self->_run_cmd("find $output_dir -name '*$pre$basename$suff'") );
        map { chomp $_; } @files;
        if ( $self->{verbose} ) { print STDERR join( "\n", @files ); print STDERR "\n"; }
        if ( $self->{verbose} ) { print STDERR "======== &runBWA: Finished ========\n"; }
        return \@files;
    }
}

=head2 bwaPostProcess

 Title   : bwaPostProcess
 Usage   : $lgtseek->bwaPostProcess(({'donor_bams' => \@donors,'host_bams' => \@hosts})
 Function: Classify the results of a short read mapping (UM, UU, MM, UM_UM etc.)
 Returns : An object with counts of the different classes as well as the path to bam files 
           containing these reads.
           $object->{files}->{           }
                            'lgt_donor'   => "$self->{output_dir}/".$prefix."lgt_donor.bam",
                            'lgt_host'    => "$self->{output_dir}/".$prefix."lgt_host.bam",
                            'integration_site_donor_donor' => "$self->{output_dir}/".$prefix."integration_site_donor_donor.bam",
                            'integration_site_donor_host' => "$self->{output_dir}/".$prefix."integration_site_donor_host.bam",
                            'microbiome_donor' => "$self->{output_dir}/".$prefix."microbiome.bam",
            $object->{counts}->{lgt}
            $object->{counts}->{microbiome}

 Args    : An object with donor and optionally host bam files.

=cut

sub bwaPostProcess {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &bwaPostProcess: Start ========\n"; }
    my $retval;

    # Do we have both donor and host bams?
    if ( $config->{donor_bams} && $config->{host_bams} ) {
        $retval = $self->_bwaPostProcessDonorHostPaired($config);
    }
    if ( $self->{verbose} ) { print STDERR "======== &bwaPostProcess: Finished ========\n"; }
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
    my ( $self, $config ) = @_;

    my @donor_fh;
    my @donor_head;

    $self->{samtools_bin} = $self->{samtools_bin} ? $self->{samtools_bin} : 'samtools';
    my $samtools = $self->{samtools_bin};
    my $prefix = $config->{output_prefix} ? "$config->{output_prefix}_" : '';

    # Open all the donor files
    map {
        if ( $self->{verbose} ) { print STDERR "Opening $_\n"; }
        if ( $_ =~ /.bam$/ ) {
            push( @donor_head, `$samtools view -H $_` );
            open( my $fh, "-|", "$samtools view $_" );
            push( @donor_fh, $fh );
        }
        elsif ( $_ =~ /.sam.gz$/ ) {
            push( @donor_head, `zcat $_ | $samtools view -H -S -` );
            open( my $fh, "-|", "zcat $_ | $samtools view -S -" );
            push( @donor_fh, $fh );
        }
    } @{ $config->{donor_bams} };

    my $class_to_file_name = {
        'lgt_donor'        => "$self->{output_dir}/" . $prefix . "lgt_donor.bam",
        'microbiome_donor' => "$self->{output_dir}/" . $prefix . "microbiome.bam"
    };

    # Check if these files exist already. If they do we'll skip regenerating them.
    my $files_exist = 1;

    map {
        if ( !-e $class_to_file_name->{$_} ) {
            $files_exist = 0;
        }
    } keys %$class_to_file_name;

    my $class_counts = {
        'lgt'        => undef,
        'microbiome' => undef
    };

    # If the files are already there and we aren't being forced to overwrite, we'll
    # just get the counts and return
    if ( $files_exist && !$config->{overwrite} ) {

        map {
            $_ =~ /^.*(.*)\_\w+$/;    ## KBS Added "." after "^" to stop "matches null string many times in regex"
            my $class = $1;
            if ( !$class_counts->{$class} ) {
                my $count = $self->_run_cmd("$samtools view $class_to_file_name->{$_} | wc -l");
                chomp $count;
                $class_counts->{$class} = $count;
                if ( $self->{verbose} ) { print STDERR "$count for $class\n"; }
            }
        } keys %$class_to_file_name;
        return {
            counts => $class_counts,
            files  => $class_to_file_name
        };
    }

    # Here are a bunch of file handles we'll use later.
    open( my $lgtd, "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "lgt_donor.bam -" )
        or die "Unable to open\n";
    open( my $microbiome_donor, "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "microbiome.bam -" )
        or die "Unable to open\n";

    my $class_to_file = {

    };

    my $more_lines = 1;

    while ($more_lines) {

        my @donor_lines;

        my $dr1_line;
        my $dr2_line;

        my $obj = $self->_getPairedClass( { fhs => \@donor_fh } );
        my $class = $obj->{class};
        $more_lines = $obj->{more_lines};
        $dr1_line   = $obj->{r1_line};
        $dr2_line   = $obj->{r2_line};
    }
}

sub _getPairedClass {
    my ( $self, $config ) = @_;

    my $fhs        = $config->{fhs};
    my $more_lines = 1;

    my $r1_class;
    my $r1_line;
    my $r2_class;
    my $r2_line;

    # Next establish the class of the donor read
    foreach my $fh (@$fhs) {
        my $r1 = <$fh>;
        my $r2 = <$fh>;
        if ( $config->{strip_xa} ) {
            chomp $r1;
            $r1 =~ s/\tXA:Z\S+$//;
            chomp $r2;
            $r2 =~ s/\tXA:Z\S+$//;
        }

        #         print STDERR "Processing $hr1$hr2$dr1$dr2";
        # Should check if these ended at the same time?
        if ( !$r1 || !$r2 ) {
            $more_lines = 0;
            last;
        }

        my $r1_flag = $self->_parseFlag( ( split( /\t/, $r1 ) )[1] );

        #            my $dr2_flag = $self->_parseFlag((split(/\t/,$dr2))[1]);
        if ( !$r1_flag->{'qunmapped'} ) {
            $r1_line  = $r1;
            $r1_class = 'M';
        }
        elsif ( !$r1_class ) {
            $r1_line  = $r1;
            $r1_class = 'U';
        }
        if ( !$r1_flag->{'munmapped'} ) {
            $r2_line  = $r2;
            $r2_class = 'M';
        }
        elsif ( !$r2_class ) {
            $r2_line  = $r2;
            $r2_class = 'U';
        }
    }

    my $class = "$r1_class$r2_class";

    return {
        class      => $class,
        r1_line    => $r1_line,
        r2_line    => $r2_line,
        more_lines => $more_lines
    };
}

sub _bwaPostProcessDonorHostPaired {
    my ( $self, $config ) = @_;

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
        'lgt_donor'                    => "$self->{output_dir}/" . $prefix . "lgt_donor.bam",
        'lgt_host'                     => "$self->{output_dir}/" . $prefix . "lgt_host.bam",
        'integration_site_donor_donor' => "$self->{output_dir}/" . $prefix . "integration_site_donor_donor.bam",
        'integration_site_donor_host'  => "$self->{output_dir}/" . $prefix . "integration_site_donor_host.bam",
        'microbiome_donor'             => "$self->{output_dir}/" . $prefix . "microbiome.bam",
    };

    # Check if these files exist already. If they do we'll skip regenerating them.
    my $files_exist = 1;

    map {
        if ( !-e $class_to_file_name->{$_} ) {
            $files_exist = 0;
        }
    } keys %$class_to_file_name;

    my $class_counts = {
        'lgt'                    => undef,
        'integration_site_host'  => undef,
        'integration_site_donor' => undef,
        'microbiome'             => undef,
        'host'                   => undef,
        'no_map'                 => undef,
        'all_map'                => undef,
        'single_map'             => undef
    };

    # If the files are already there and we aren't being forced to overwrite, we'll
    # just get the counts and return
    if ( $files_exist && !$config->{overwrite} ) {

        map {
            $_ =~ /^.*(.*)\_\w+$/;    ## KBS added "." after ^ to stop "matches null string many times in regex"
            my $class = $1;
            if ( !$class_counts->{$class} ) {
                my $count = $self->_run_cmd("$samtools view $class_to_file_name->{$_} | cut -f1 | uniq | wc -l");
                chomp $count;
                $class_counts->{$class} = $count * 1;
                if ( $self->{verbose} ) { print STDERR "$count for $class\n"; }
            }
        } keys %$class_to_file_name;
        my $total = $self->_run_cmd( "$samtools view " . $config->{donor_bams}[0] . " | cut -f1 | uniq | wc -l" );
        chomp $total;
        $class_counts->{total} = $total * 1;
        return {
            counts => $class_counts,
            files  => $class_to_file_name
        };
    }

    # Here are a bunch of file handles we'll use later.
    if ( $self->{verbose} ) { print STDERR "$self->{output_dir}/" . $prefix . "lgt_donor.bam\n"; }
    open( my $lgtd, "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "lgt_donor.bam -" )
        or die "Unable to open\n";
    open( my $lgth, "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "lgt_host.bam -" )
        or die "Unable to open\n";
    open( my $int_site_donor_d,
        "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "integration_site_donor_donor.bam -" )
        or die "Unable to open\n";
    open( my $int_site_donor_h,
        "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "integration_site_donor_host.bam -" )
        or die "Unable to open\n";
    open( my $microbiome_donor, "| $samtools view -S -b -o $self->{output_dir}/" . $prefix . "microbiome.bam -" )
        or die "Unable to open\n";

    my $class_to_file = {
        'lgt_donor'                    => $lgtd,
        'lgt_host'                     => $lgth,
        'integration_site_donor_donor' => $int_site_donor_d,
        'integration_site_donor_host'  => $int_site_donor_h,
        'microbiome_donor'             => $microbiome_donor
    };

    my @donor_fh;
    my @host_fh;
    my @donor_head;
    my @host_head;

    # Open all the donor files
    map {
        if ( $self->{verbose} ) { print STDERR "Opening $_\n"; }

        if ( $_ =~ /.bam$/ ) {
            push( @donor_head, `$samtools view -H $_` );
            open( my $fh, "-|", "$samtools view $_" );
            push( @donor_fh, $fh );
        }
        elsif ( $_ =~ /.sam.gz$/ ) {
            push( @donor_head, `zcat $_ | $samtools view -H -S -` );
            open( my $fh, "-|", "zcat $_ | $samtools view -S -" );
            push( @donor_fh, $fh );
        }
    } @{ $config->{donor_bams} };

    # Open all the host files
    map {
        if ( $self->{verbose} ) { print STDERR "Opening $_\n"; }
        if ( $_ =~ /.bam$/ ) {
            push( @host_head, `$samtools view -H $_` );
            open( my $fh, "-|", "$samtools view $_" );
            push( @host_fh, $fh );
        }
        elsif ( $_ =~ /.sam.gz$/ ) {
            push( @host_head, `zcat $_ | $samtools view -H -S -` );
            open( my $fh, "-|", "zcat $_ | $samtools view -S -" );
            push( @host_fh, $fh );
        }
    } @{ $config->{host_bams} };

    # Prime the files with headers.
    map {
        if ( $_ =~ /_donor$/ ) {
            print "Printing header to $_ donor file\n";
            print { $class_to_file->{$_} } join( '', @donor_head );

            #            close $class_to_file->{$_};
        }
        elsif ( $_ =~ /_host$/ ) {
            print "Printing header to $_ host file\n";
            print { $class_to_file->{$_} } join( '', @host_head );

            #            close $class_to_file->{$_};
        }
    } keys %$class_to_file;

    #   exit;
    my $more_lines = 1;

    my $line_num = 0;

    while ($more_lines) {

        my @donor_lines;
        my @host_lines;

        # Get the class of the host mappings
        my $obj = $self->_getPairedClass( { fhs => \@host_fh } );
        my $hclass = $obj->{class};
        $more_lines = $obj->{more_lines};
        my $hr1_line = $obj->{r1_line};
        my $hr2_line = $obj->{r2_line};

        # Get the class of the donor mappings
        $obj = $self->_getPairedClass( { fhs => \@donor_fh } );
        my $dclass = $obj->{class};
        $more_lines = $obj->{more_lines};
        my $dr1_line = $obj->{r1_line};
        my $dr2_line = $obj->{r2_line};

        if ($more_lines) {
            my $paired_class = "$dclass\_$hclass";

            # print the donor lines to the donor file (if we are keeping this output file)
            if ( $class_to_file->{ $classes_both->{$paired_class} . "_donor" } ) {
                print { $class_to_file->{ $classes_both->{$paired_class} . "_donor" } } "$dr1_line\n$dr2_line\n";
            }

            # print the host lines to the host file (if we are keeping this output file)
            if ( $class_to_file->{ $classes_both->{$paired_class} . "_host" } ) {
                print { $class_to_file->{ $classes_both->{$paired_class} . "_host" } } "$hr1_line\n$hr2_line\n";
            }

            # Increment the count for this class
            if ( $classes_both->{$paired_class} ) {
                $class_counts->{ $classes_both->{$paired_class} }++;
            }

            # Increment the total count
            $line_num++;
        }
    }

    # Close up the file handles
    map {
        if ( $_ =~ /_donor$/ ) {
            if ( $self->{verbose} ) { print STDERR "closing $_ donor file\n"; }
            close $class_to_file->{$_};
        }
        elsif ( $_ =~ /_host$/ ) {
            if ( $self->{verbose} ) { print STDERR "closing $_ host file\n"; }
            close $class_to_file->{$_};
        }
    } keys %$class_to_file;

    # Set the total
    $class_counts->{total} = $line_num;

    # Return the files and the counts
    return {
        counts => $class_counts,
        files  => $class_to_file_name
    };
}

=head2 blast2lca

 Title   : blast2lca
 Usage   : my $lca_file = $lgtseek->blast2lca({'blast' => 'path_to_blast-m8.txt','output_dir' => 'path_to_out_dir'})
 Function: Take a blast -m8 report and calculate the LCA'
 Returns : A object with the path to the output file(s) with LCA's.                             
            $lca_file->{independent} = Path to the file with the LCA's for each unique ID. 
            $lca_file->{PE_lca}      = Path to the file with the conservative and liberal LCA for combining mate information
 Args    : 
            blast               => 
            output_dir          =>
            evalue_cutoff       =>  ###  [1] Max evalue allowed for a hit. Example : 1e-5.
            best_hits_only      => <0|1> [0] 1= Parse the Blast file for only best hits.
            combine_PE_lca      => <0|1> [0] 1= Merge the PE reads with the conservative and liberal methods. ## NOT IMPLEMENTED YET ##

=cut

sub blast2lca {
    my ( $self, $config ) = @_;

    my $gi2tax = $self->getGiTaxon( {} );
    if ( $self->{verbose} ) { print STDERR "======== &blast2lca: Start ========\n"; }
    open IN, "<$config->{blast}" or confess "Unable to open $config->{blast}\n";
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");

    my ( $name, $directories, $suffix ) = fileparse( $config->{blast}, qr/\.[^.]*/ );
    $name = $config->{output_prefix} ? $config->{output_prefix} : $name;
    open OUT, ">$output_dir/$name\_lca-independent.out"
        or confess "Error: &blast2lca unable to open $output_dir/$name\_lca-independent.out";
    if ( $config->{best_hits_only} == 1 ) {
        print OUT join( "\t", ( 'read_id', 'lca', 'best_evalue' ) );
    }
    else {
        print OUT join( "\t", ( 'read_id', 'lca', 'best_evalue', 'highest_evalue' ) );
    }    ## KBS
    print OUT "\n";
    my $hits_by_readname = {};
    my $id;
    my $evalue_cutoff = 1;
    my $best_evalue;
    my $worst_evalue;
    my $lca;

    my @vals;
    while (<IN>) {
        my @fields = split(/\t/);
        my $new_id = $fields[0];
        if ( $fields[0] =~ /^(.*)(\/\d)*/ ) {
            $new_id = $1;
        }
        if ( $config->{evalue_cutoff} ) { $evalue_cutoff = $config->{evalue_cutoff}; }
        next if ( $new_id eq $id && $fields[10] > $evalue_cutoff );
        my $taxon = $gi2tax->getTaxon( $fields[1] );
        $taxon->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;
        if ( !defined($id) && $fields[10] < $evalue_cutoff ) {
            $id = $new_id;
            if ( $config->{best_hits_only} == 1 ) { $evalue_cutoff = $fields[10]; }
            $best_evalue                                              = $fields[10];
            $worst_evalue                                             = $fields[10];
            $lca                                                      = $taxon->{lineage};
            $hits_by_readname->{ $fields[0] }->{ $taxon->{taxon_id} } = $taxon->{lineage};
        }
        elsif ( $new_id eq $id && $fields[10] <= $evalue_cutoff ) {
            $hits_by_readname->{ $fields[0] }->{ $taxon->{taxon_id} } = $taxon->{lineage};
            my $newlca = &_find_lca( [ $lca, $taxon->{lineage} ] );
            $lca = $newlca;
            if ( $fields[10] > $worst_evalue ) { $worst_evalue = $fields[10]; }
        }
        elsif ( $new_id ne $id ) {
            if ( $config->{best_hits_only} == 1 ) {
                print OUT "$id\t$lca\t$best_evalue\n";
            }
            else {
                print OUT "$id\t$lca\t$best_evalue\t$worst_evalue\n";
            }
            $id = $new_id;
            if ( $config->{best_hits_only} == 1 ) { $evalue_cutoff = $fields[10]; }
            $best_evalue  = $fields[10];
            $worst_evalue = $fields[10];
            $lca          = $taxon->{lineage};
        }
    }

    # Print out last line
    if ( $config->{best_hits_only} == 1 ) {
        print OUT "$id\t$lca\t$best_evalue\n";
    }
    else {
        print OUT "$id\t$lca\t$best_evalue\t$worst_evalue\n";
    }
    close OUT;
    if ( $config->{combine_PE_lca} == 1 ) {
        open( IN, "<", "$output_dir/$name\_lca-independent.out" )
            or confess "Error: &blast2lca unable to open input: $output_dir/$name\_lca-independent.out because: $!\n";
        open( OUT1, ">", "$output_dir/$name\_lca-conservative.out" )
            or confess "Error: &blast2lca unable to open output: $output_dir/$name\_lca-conservative.out because: $!\n";
        open( OUT2, ">", "$output_dir/$name\_lca-liberal.out" )
            or confess "Error: &blast2lca unable to open output: $output_dir/$name\_lca-liberal.out because: $!\n";
        my $continue = 1;
        while ( $continue == 1 ) {
            my $line1 = <IN>;
            my $line2 = <IN>;
            if ( !$line1 || !$line2 ) { $continue = 0; }
        NEXT:
            my @f1 = split( /\t/, $line1 );
            my @f2 = split( /\t/, $line2 );
            my $id1 = $f1[0];
            my $id2 = $f2[0];
            $id1 =~ /([A-Za-z0-9-.|:]+)(\_?[1,2]?)/;
            my $id1_short = $1;
            $id2 =~ /([A-Za-z0-9-.|:]+)(\_?[1,2]?)/;
            my $id2_short = $1;

            if ( $id1_short ne $id2_short ) {
                $line1 = $line2;
                $line2 = <IN>;
                goto NEXT;
            }
            my $lca1             = @f1[1];
            my $lca2             = @f2[1];
            my $conservative_lca = &_find_lca( [ $lca1, $lca2 ] );

            # Calculate Liberal LCA
            my $liberal_lca;
            if ( $lca1 =~ $lca2 || $lca2 =~ $lca1 ) {
                if ( length($lca1) >= length($lca2) ) {
                    $liberal_lca = $lca1;
                }
                else {
                    $liberal_lca = $lca2;
                }
            }
            print OUT1 "$id1_short\t$conservative_lca\n";
            print OUT2 "$id1_short\t$liberal_lca\n";
        }
        close OUT1;
        close OUT2;
    }
    if ( $self->{verbose} ) { print STDERR "======== &blast2lca: Finished ========\n"; }
    return "$output_dir/$name\_lca-independent.out";
}

sub _prep_lca_print {
    my $hits_by_readname = shift;
    my $vals             = shift;

    # Find the LCA of only mated hits
    my @reads = keys %$hits_by_readname;
    if ( scalar @reads > 1 ) {

        # First we'll find the LCA requiring both reads to match the same taxon id.
        my @good_lineages;
        my @read1_lineages;
        map {
            if ( $hits_by_readname->{ $reads[1] }->{$_} ) {
                push( @good_lineages, $hits_by_readname->{ $reads[1] }->{$_} );
            }
            push( @read1_lineages, $hits_by_readname->{ $reads[0] }->{$_} );
        } keys %{ $hits_by_readname->{ $reads[0] } };

        my @read2_lineages;
        map { push( @read2_lineages, $hits_by_readname->{ $reads[1] }->{$_} ); }
            keys %{ $hits_by_readname->{ $reads[1] } };

        my $paired_lca = &_find_lca( \@good_lineages );
        push( @$vals, $paired_lca );

        # Next we'll find the LCA where if either LCA is a substring of the
        # other we'll take the longer one (more specific).
        my $read1_lca = &_find_lca( \@read1_lineages );
        my $read2_lca = &_find_lca( \@read2_lineages );

        ## These two LCA's only apply to combine PE read LCA's
        my $conservative_lca = &_find_lca( [ $read1_lca, $read2_lca ] );  ## Take the shorter but most parsimonious LCA.

        my $liberal_lca;    ## Take the longer LCA if the shorter LCA is a substring of the longer LCA.
        if ( $read1_lca =~ $read2_lca || $read2_lca =~ $read1_lca ) {
            if ( length($read1_lca) >= length($read2_lca) ) { $liberal_lca = $read1_lca; }
            if ( length($read2_lca) >= length($read1_lca) ) { $liberal_lca = $read2_lca; }
        }

        push( @$vals, $liberal_lca );
        push( @$vals, $conservative_lca );
        $hits_by_readname = {};
    }
    else {
        push( @$vals, ( '', '' ) );
    }
    return $vals;
}

sub _find_lca {
    my $lineages = shift;

    # prime LCA
    my @lca = split( ';', $lineages->[0] );

    foreach my $l (@$lineages) {
        my $newlca = [];
        my @lineage = split( ';', $l );
        for ( my $i = 0; $i < @lineage; $i++ ) {
            if ( $lca[$i] eq $lineage[$i] ) {
                push( @$newlca, $lineage[$i] );
            }
            else {
                last;
            }
        }
        @lca = @$newlca;
    }
    return join( ';', @lca );
}

=head2 &bestBlast2

 Title   : bestBlast2
 Usage   : $lgtseek->bestBlast2(({'fasta' => 'fasta_file', 'ref' => '/path/to/ref/db'}))
 Function: Run blast (or megablast) against a reference database to find best blast hits.
 Returns : Hash with:
            overall_blast => $overallfile,
            out1file_file => $out1file,
            out2file_file => $out2file
 Args    : An object with the input bame files and the reference database.
            fasta =>
            bam => 
            lineage1 =>
            lineage2 => 
            db =>
            output_dir =>

=cut

sub bestBlast2 {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) { print STDERR "======== &bestBlast2: Start ========\n"; }
    my $output_dir = defined $config->{output_dir} ? "$config->{output_dir}" : "$self->{output_dir}";
    $self->_run_cmd("mkdir -p $output_dir");

    # Convert bams to fasta
    my $fasta;
    if ( $config->{'bam'} ) {
        my $newfasta = $self->sam2Fasta(
            {   input      => $config->{'bam'},
                output_dir => $config->{output_dir},
                paired     => 1
            }
        );
        $fasta = $newfasta;
    }
    elsif ( $config->{'fasta'} ) {
        $fasta = $config->{fasta};
    }

    # Blast fasta @ database and filter for best hits.
    my $files = LGTBestBlast::filterBlast(
        {   blast_bin  => $config->{blast_bin},
            fasta      => $fasta,
            db         => $config->{db},
            gitaxon    => $self->getGiTaxon( {} ),
            lineage1   => $config->{lineage1},
            lineage2   => $config->{lineage2},
            output_dir => $output_dir,
        }
    );

# my ($name,$directories,$suffix) = fileparse($fasta,qr/\.[^.]*/);
# open OUT, ">$self->{output_dir}/$name\_filtered_blast.list" or die "Unable to open $self->{output_dir}/$name\_filtered_blast.list\n";
# my @outputs;
# map {
# print OUT "$files->{$_}\n";
# } keys %$files;

    # close OUT;

    # $files->{list_file} = "$self->{output_dir}/$name\_filtered_blast.list";
    if ( $self->{verbose} ) { print STDERR "======== &bestBlast2: Finished ========\n"; }
    return $files;
}

sub OUTDATED_bestBlast {
    my ( $self, $config ) = @_;

    my @fastas;

    # Convert bams to fasta
    if ( $config->{'bams'} ) {

    }
    $self->{output_dir} = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    if ( $config->{'fastas'} ) {
        push( @fastas, @{ $config->{'fastas'} } );
    }
    my $blast_bin = '';
    if ( $config->{blast_bin} ) {
        $blast_bin = $config->{"blast_bin"} . " -d $config->{'blast_db'} -e 1 -m8 -T F ";
    }
    my @overall_files;
    my @lineage1_files;
    my @lineage2_files;
    foreach my $fasta (@fastas) {

        my $cmd
            = "perl -I $self->{bin_dir} $self->{bin_dir}/filter_lgt_best_hit.pl"
            . " --input=$fasta"
            . " --blast_bin=\'$blast_bin\'"
            . " --taxon_dir="
            . $self->{taxon_dir}
            . " --trace_mapping="
            . $config->{trace_mapping}
            . " --lineage1="
            . $config->{lineage1}
            . " --lineage2="
            . $config->{lineage2}
            . " --filter_lineage="
            . $config->{filter_lineage}
            . " --dbhost="
            . $self->{taxon_host}
            . " --idx_dir="
            . $self->{taxon_idx_dir}
            . " --output_dir="
            . $self->{output_dir};

        $self->_run_cmd($cmd);
        my ( $name, $directories, $suffix ) = fileparse( $fasta, qr/\.[^.]*/ );

        push( @overall_files,  "$self->{output_dir}/$name\_overall.out" );
        push( @lineage1_files, "$self->{output_dir}/$name\_lineage1.out" );
        push( @lineage2_files, "$self->{output_dir}/$name\_lineage2.out" );
    }

    # Should merge the three file types here.
    my $cmd = "cat " . join( ' ', @overall_files ) . " > $self->{output_dir}/all_overall.out";
    $self->_run_cmd($cmd);
    $cmd = "cat " . join( ' ', @lineage1_files ) . " > $self->{output_dir}/all_lineage1.out";
    $self->_run_cmd($cmd);
    $cmd = "cat " . join( ' ', @lineage2_files ) . " > $self->{output_dir}/all_lineage2.out";
    $self->_run_cmd($cmd);

    open OUT, ">$self->{output_dir}/filtered_blast.list"
        or die "Unable to open $self->{output_dir}/filtered_blast.list\n";
    print OUT
        "$self->{output_dir}/all_overall.out\n$self->{output_dir}/all_lineage1.out\n$self->{output_dir}/all_lineage2.out";
    close OUT;

    return {
        overall_blast  => "$self->{output_dir}/all_overall.out",
        lineage1_blast => "$self->{output_dir}/all_lineage1.out",
        lineage2_blast => "$self->{output_dir}/all_lineage2.out",
        list_file      => "$self->{output_dir}/filtered_blast.list"
    };
}

=head2 &runLgtFinder

 Title   : runLgtFinder
 Usage   : $lgtseek->runLgtFinder(({'inputs' => \@files})
 Function: Take output from LGTBestBlast and look for putative LGT's both within reads and across read pairs
 Returns : A file with all the LGT's hit information
 Args    : An object from BestBlast2
            max_overlap => 20,
            min_length  => 0,
            ref_lineage => 'Homo',
            output_dir  => '/somewhere/'

=cut

sub runLgtFinder {
    my ( $self, $config ) = @_;
    if ( $config->{output_dir} ) {
        $self->_run_cmd("mkdir -p $config->{output_dir}");
    }
    LGTFinder::findLGT($config);
}

=head2 splitBam

 Title   : splitBam
 Usage   : $lgtseek->splitBam(({'input' => $file})
 Function: Split a bam file into smaller chunks.
 Returns : A list of the bam files split up
 Args    : 
            input = An object with the input bam files
            seqs_per_file = # of seqs per file for each split file
            output_dir = Directory for output
            samtools_bin = bin directory with samtools

=cut

sub splitBam {
    my ( $self, $config ) = @_;
    if ( !$config->{input} ) { die "Must give &splitBam an input =>\n"; }
    if ( $self->empty_chk( { input => $config->{input} } ) == 1 ) {
        print STDERR "Warning: &split_bam input: $config->{input} is empty.\n";
    }
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) { print STDERR "======== &splitBam: Start ========\n"; }
    $self->_run_cmd("mkdir -p $output_dir");
    my $seqs_per_file = $config->{seqs_per_file} ? $config->{seqs_per_file} : 10000;
    my $samtools = $self->{samtools_bin};

    # Parse out the pieces of the filename
    my ( $fn, $path, $suffix ) = fileparse( $config->{input}, '.bam' );

    # Pull out the header
    my $header = $self->_run_cmd("$samtools view -H $config->{input}");

    # Open the input file
    open( IN, "-|", "$samtools view $config->{input}" );

    # Open the first output file
    my $count = 0;
    my $ofile = "$output_dir/$fn\_$count.bam";
    open( my $ofh, "| samtools view -S -b -o $ofile -" );
    my @outfiles = ($ofile);

    # Print the header
    print $ofh $header;

    my $i = 0;
    while ( my $line = <IN> ) {
        my @fields = split( /\t/, $line );
        my $flag = $self->_parseFlag( $fields[1] );

        # Strip out the XA tag
        $line =~ s/\s+XA:Z:\S+//;

        # Make sure we don't accidentally split read pairs.
        if ( $flag->{'last'} && $i >= $seqs_per_file ) {

            # Print out the line
            print $ofh $line;

            # Close the old file
            close $ofh;

            # Open the new file
            $count++;
            $ofile = "$output_dir/$fn\_$count.bam";
            open( $ofh, "| samtools view -S -b -o $ofile -" );
            push( @outfiles, $ofile );

            # Print the header to the new file
            print $ofh $header;

            # Reset the counter
            $i = 0;
        }

        # If we haven't filled the current file yet just keep printing.
        print $ofh $line;

        # Increment the counter
        $i++;
    }
    close $ofh;
    if ( $self->{verbose} ) {
        print STDERR "Split $config->{input} into $count bams, each with $seqs_per_file sequences per bam.\n";
    }
    my $files = \@outfiles;
    if ( $self->{verbose} ) { print STDERR "======== &splitBam: Finished ========\n"; }
    return $files;
}

sub _dec2bin {
    my $self = shift;
    my $str = unpack( "B32", pack( "N", shift ) );
    $str =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
    return $str;
}

sub _parseFlag {
    my $self   = shift;
    my $int    = shift;
    my $rawbin = $self->_dec2bin($int);
    my $rev    = scalar $rawbin;
    if ( $rev eq $rawbin ) {

        #    print "ERROR $rev $rawbin\n";
    }
    my $bin = sprintf( "%011d", $rev );
    my $final_bin = reverse $bin;
    return {
        'paired'    => substr( $final_bin, 0,  1 ),
        'proper'    => substr( $final_bin, 1,  1 ),
        'qunmapped' => substr( $final_bin, 2,  1 ),
        'munmapped' => substr( $final_bin, 3,  1 ),
        'qrev'      => substr( $final_bin, 4,  1 ),
        'mrev'      => substr( $final_bin, 5,  1 ),
        'first'     => substr( $final_bin, 6,  1 ),
        'last'      => substr( $final_bin, 7,  1 ),
        'secondary' => substr( $final_bin, 8,  1 ),
        'failqual'  => substr( $final_bin, 9,  1 ),
        'pcrdup'    => substr( $final_bin, 10, 1 )
    };
}

=head2 decrypt

 Title   : decrypt
 Usage   : my $decrypted_bam=$lgtseek->decrypt(({'input' => <file.bam.gpg>, 'key' => <Path to key>})
 Function: Decrypt a .bam.gpg input file. 
 Returns : An un-encrypted bam
 Args    : 
    input       => encrypted bam for decryption
    url         => url to download the key from
    key         => path to key file
    output_dir  => directory for output

=cut

sub decrypt {
    my ( $self, $config ) = @_;
    if ( !$config->{input} || $config->{input} !~ /\.bam\.gpg$/ ) {
        die "Error: Must pass an encrypted .bam.gpg to this subroutine.\n";
    }
    if ( !$config->{url} && !$config->{key} ) {
        die "Must give an url to download the key from or the path to the key.\n";
    }
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) { print STDERR "======== &Decrypt: Start ========\n"; }
    $self->_run_cmd("mkdir -p $output_dir");
    my ( $bam, $path, $suffix ) = fileparse( $config->{input}, '.gpg' );
    my $key;
    if ( $config->{url} ) {
        $self->_run_cmd("wget $config->{url}");
        $key = `find . -name '*.key' | cut -f2 -d '/'`;
    }
    else {
        $key = $config->{key};
    }
    $self->_run_cmd("gpg --import $key");
    $self->_run_cmd("gpg -o $output_dir/$bam -d $config->{input}");
    my $outfile = "$output_dir/$bam";
    $self->_run_cmd("rm $config->{input}");
    $self->_run_cmd("rm $key");
    if ( $self->{verbose} ) { print STDERR "======== &Decrypt: Finished ========\n"; }
    return $outfile;
}

=head2 mpileup

 Title   : mpileup
 Usage   : my $mpileup_file=$lgtseek->mpileup({'input' => <bam>})
 Function: Calculate samtools mpileup 
 Returns : File path to mpileup output.
 Args    : 
    input       => unsorted bam
    srtd_bam     => position sorted bam input
    ref         =>Reference (Optional)
    max         => max per-BAM depth to avoid excessive memory usage [250] (Optional)
    cleanup     => <0|1> [0] 1= rm input.srt.bam 
    overwrite   => <0|1> [0] 1= overwrite 

=cut

sub mpileup {
    my ( $self, $config ) = @_;
    if ( !$config->{input} && !$config->{srtd_bam} ) { die "Error: Must give an input bam to calculat coverage on.\n"; }
    my $input = $config->{input} ? $config->{input} : $config->{srtd_bam};
    if ( $self->empty_chk( { input => $input } ) == 1 ) { print STDERR "Warning: &mpileup input: $input is empty.\n"; }
    my $overwrite = $config->{overwrite} ? $config->{overwrite} : 0;
    my ( $fn, $path, $suffix ) = fileparse( $config->{input}, '.bam' );
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    my $output = "$output_dir/$fn\.mpileup";
    if ( $self->{verbose} ) { print STDERR "======== &mpileup: Start ========\n"; }

    if ( -e $output && $overwrite == 0 ) {
        print STDERR "Already found the output for &mpileup: $output\n";
        return $output;
    }
    $self->_run_cmd("mkdir -p $output_dir");
    my $samtools = $self->{samtools_bin};
    my $srtd_bam;
    my $max_opt = $config->{max}     ? "-d $config->{max}" : "";
    my $ref_opt = $config->{ref}     ? "-f $config->{ref}" : "";
    my $cleanup = $config->{cleanup} ? $config->{cleanup}  : "0";

    my $cmd = "samtools";

    ## New Method
    if ( $config->{input} ) {
        $cmd      = $cmd . "sort -m 5000000000 -o $config->{input} - | samtools mpileup";
        $srtd_bam = "-";
    }
    else {
        $cmd      = $cmd . "mpileup";
        $srtd_bam = $config->{srtd_bam};
    }
    $cmd = $cmd . "-A $max_opt $ref_opt $srtd_bam > $output";
    $self->_run_cmd($cmd);
    if ( $self->{verbose} ) { print STDERR "======== &mpileup: Finished ========\n"; }
    return $output;
}

=head2 &prelim_filter   

 Title   : prelim_filter
 Usage   : my $potential_LGT_bam = $lgtseek->prelim_filter({input => <bam>})
 Function: Removes M_M reads
 Returns : A bam will all other reads
 Args    : A hash containing potentially several config options:
        input_bam       => bam
        output_dir      => directory for output
        keep_softclip   => <0|1> [0] 1= Keep the soft clipped M_M reads 
        overwrite       => <0|1> [1] 1= Overwrite the output if it already exists
        split_bam       => <0|1> [0] 1= Split by by #seqs_per_file
        seqs_per_file   => [50000000] 
        name_sort_input => <0|1> [0] 1= Resort input bam on names

=cut

sub prelim_filter {

    ## General code scheme: (1) Resort based on names. (2) Filter out M_M reads. (3) Then finally split into chunks.
    my ( $self, $config ) = @_;
    my $input = $config->{input_bam} ? $config->{input_bam} : $self->{input_bam};
    if ( $self->empty_chk( { input => $input } ) == 1 ) {
        print STDERR "Warning: &prelim_filter input: $input is empty.\n";
    }
    my $output_dir = $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");

    my $keep_softclip   = defined $config->{keep_softclip}   ? $config->{keep_softclip}   : $self->{keep_softclip};
    my $overwrite       = defined $config->{overwrite}       ? $config->{overwrite}       : $self->{overwrite};
    my $split_bam       = defined $config->{split_bam}       ? $config->{split_bam}       : $self->{split_bam};
    my $prelim_filter   = defined $config->{prelim_filter}   ? $config->{prelim_filter}   : $self->{prelim_filter};
    my $seqs_per_file   = $config->{seqs_per_file}           ? $config->{seqs_per_file}   : $self->{seqs_per_file};
    my $name_sort_input = defined $config->{name_sort_input} ? $config->{name_sort_input} : $self->{name_sort_input};
    if ( $self->{verbose} ) { print STDERR "======== &prelim_filter: Start ========\n"; }
    my ( $fn, $path, $suffix ) = fileparse( $input, ( '.srt.bam', '.bam' ) );
    my $header = $self->_run_cmd("samtools view -H $input");
    my @output_list;    ## Array of output bams to return
    ## Check if the output exists already
    ## $output isn't set right at this point in the scrip to check if the final output already exist. FIX later.
    my $files_exist = 0;
    if ( -e "$output_dir$fn\_0_prelim.bam" ) { $files_exist = 1; }
    if ( $files_exist == 1 && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR "Already found the output for &prelim_filter: $output_dir$fn\_0_prelim.bam\n";
        }
        chomp( my @output_list = `find $output_dir -name '$fn\*_prelim.bam'` );
        return \@output_list;
    }
    ## Check we dont' have the name-sorted.bam already
    my $sorted_bam = "$output_dir$fn\_name-sorted.bam";
    if ( -e "$sorted_bam" ) { $files_exist = 1; }
    if ( $files_exist == 1 && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR "Already found &prelim_filter \$sorted_bam, starting to filter: $sorted_bam.\n";
        }
        $prelim_filter = 0;
    }
    ## (1). Sort the bam by name instead of position
    if ( $name_sort_input == 1 ) {
        if ( $self->{verbose} ) { print STDERR "======== &prelim_filter Sort ========\n"; }
        my $cmd = "samtools sort -n -@ $self->{threads} -m $self->{sort_mem} $input $output_dir$fn\_name-sorted";
        if ( $self->{verbose} ) { print STDERR "======== &prelim_filter Sort: $cmd ===\n"; }
        $self->_run_cmd("$cmd");
    }
    else {
        $sorted_bam = $input;
    }
    ## (2). Prelim filtering.
    my $more_lines     = 1;
    my $num_pass       = 0;
    my $num_null       = 0;
    my $num_singletons = 0;
    if ( $prelim_filter == 1 ) {
        if ( $self->{verbose} ) { print STDERR "======== &prelim_filter: Filter ========\n"; }
        ## Setup and open output
        my $i      = 0;                                      ## Used to count number of lines per bam as being processed
        my $count  = 0;                                      ## Used to count number of split output bams
        my $output = "$output_dir/$fn\_$count\_prelim.bam";
        open( my $out, "| samtools view -S - -bo $output" ) || $self->fail("Can't open: $output because: $!\n");
        print $out "$header";
        push( @output_list, $output );

        ## Open input bam
        open( my $infh, "samtools view $sorted_bam |" ) || $self->fail("&prelim_filter can't open the \$sorted_bam: $sorted_bam because: $1\n");
        ## (2). Actual Prelim Filtering starting:
        ##      Read through and filter bam reads, splitting output as we go
        if ( $self->{verbose} ) { print STDERR "======== &prelim_filter: Sort ========\n"; }
       
        while ( $more_lines == 1 ) {
            my $read1 = <$infh>;
            my $read2 = <$infh>;
            my $print = 0;
            ## Stop if we have empty lines
            if ( !$read1 || !$read2 ) { $more_lines = 0; last; }
            ## (3). Split output: close the current output and open a new output
            if ( $i >= $seqs_per_file && $split_bam == 1 ) {
                close $out or $self->fail("Unable to close &prelim_filter output: $output\n");
                $count += 1;
                $output = "$output_dir/$fn\_$count\_prelim.bam";
                open( $out, "| samtools view -S - -bo $output" ) || $self->fail("Can't open: $output because: $!\n");
                push( @output_list, $output );
                print $out "$header";
                $i = 0;
            }
            ## Split needed fields
            my ( $id1, $flag1, $cigar1, $sequence1 ) = ( split /\t/, $read1 )[ 0, 1, 5, 9 ];
            my ( $id2, $flag2, $cigar2, $sequence2 ) = ( split /\t/, $read2 )[ 0, 1, 5, 9 ];
            ## Check sequence id's. ## Both id's must be the same && no "null" reads
            while ( $id2 =~ /null/ ) {
                $num_null++;
                $read2 = <$infh>;
                ( $id2, $flag2, $cigar2, $sequence2 ) = ( split /\t/, $read2 )[ 0, 1, 5, 9 ];
            }
            ## Both id's must be the same.
            while ( $id1 =~ /null/ || $id1 ne $id2 )
            {    ## If id1 isn't the same as id2, skip id1. Keep the pair moving down the bam until both are the same.
                if   ( $id1 =~ /null/ ) { $num_null++; }
                else                    { $num_singletons++; }
                $read1 = $read2;
                $read2 = <$infh>;
                ( $id1, $flag1, $cigar1, $sequence1 ) = ( split /\t/, $read1 )[ 0, 1, 5, 9 ];
                ( $id2, $flag2, $cigar2, $sequence2 ) = ( split /\t/, $read2 )[ 0, 1, 5, 9 ];
            }

            #next if ($id1=~/null/ || $id2=~/null/);
            ## Convert sam flag into usable data
            my $converted_flag1 = $self->_parseFlag($flag1);
            my $converted_flag2 = $self->_parseFlag($flag2);
            ## FILTER for reads that are UM (M_UM,UM_M,UM_UM)
            if ( $converted_flag1->{'qunmapped'} || $converted_flag1->{'munmapped'} ) { $print = 1; }
            ## FILTER for soft clipped reads
            if ( $keep_softclip == 1 ) {
                map {
                    if ( $_ =~ /(\d+)M(\d+)S/ && $2 >= 24 ) { $print = 1; }
                    if ( $_ =~ /(\d+)S(\d+)M/ && $1 >= 24 ) { $print = 1; }
                } ( $cigar1, $cigar2 );
            }
            if ( $print == 1 ) {
                print $out "$read1$read2";
                $i += 2;
                $num_pass++;
            }
        }
        close $out; ## or $self->fail("Can't close out: $output because: $!\n");
        ## Check to make sure the last file isn't empty.
        if($self->empty_chk({input=>$output})==1){ 
            $self->_run_cmd("rm $output"); 
             my $remove_empty_bam = pop(@output_list);
        }
        if ( $name_sort_input == 1 && $sorted_bam =~ /name-sorted.bam$/ ) { $self->_run_cmd("rm $sorted_bam"); }
    }
    else {
        push( @output_list, $sorted_bam );
    }

    my @sort_out_list = sort @output_list;
    if ( $self->{verbose} ) {
        print STDERR "======== &prelim_filter: Singletons:$num_singletons ==\n";
        print STDERR "======== &prelim_filter: Null:$num_null ==\n";
        print STDERR "======== &prelim_filter: Pass:$num_pass ==\n";
        print STDERR "======== &prelim_filter: Finished ========\n";
    }
    return \@sort_out_list;
}

=head2 &filter_bam_by_ids

 Title   : filter_bam_by_ids
 Usage   : my $filtered_bam = $lgtseek->filter_bam_by_ids({input => <bam>, good_list=> <list of desired ids>})
           output = output_dir/prefix_suffix.bam
 Function: Filter bam by ids
 Returns : A hash{bam} = new filtered bam 
           A hash{count} = # ids found
 Args    : A hash containing potentially several config options:
        input_bam     => bam for filtering
        output_dir    => directory for output
        output_prefix => Prefix for output
        output_suffix => Suffix for output
        output        => Full Path, prefix, name, and suffix for output. 
        good_list     => File path with list of desired reads
        bad_list      => File path with a list of reads to discard
=cut

sub filter_bam_by_ids {
    my ( $self, $config ) = @_;
    if ( !$config->{input_bam} ) {
        $self->fail("Error: Must pass &filter_bam_by_ids an input_bam => <BAM_TO_FILTER>\n");
    }
    if ( $self->empty_chk( { input => "$config->{input_bam}" } ) == 1 ) {
        print STDERR "Warning: Can not filter ids from an empty input bam: $config->{input_bam}\n";
        return { count => 0, file => $config->{input_bam} };
    }
    if ( !$config->{good_list} && !$config->{bad_list} ) {
        $self->fail(
            "Error: Must pass &filter_bam_by_ids a file with a list of reads to filter on. Use good_list => || bad_list => \n"
        );
    }
    if ( $self->{verbose} ) { print STDERR "======== &filter_bam_by_ids: Start ========\n"; }
    ## Setup hash of ids
    my $good_ids = {};
    my $bad_ids  = {};
    my %found_ids;
    if    ( $config->{good_list} ) { $good_ids = $self->_read_ids( { list => $config->{good_list} } ); }
    elsif ( $config->{bad_list} )  { $bad_ids  = $self->_read_ids( { list => $config->{bad_list} } ); }
    ## Setup input and output bams
    my $input = $config->{input_bam};
    my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
    my $out_dir = defined $config->{output_dir}    ? $config->{output_dir}    : $path;
    my $prefix  = defined $config->{output_prefix} ? $config->{output_prefix} : $fn;
    my $suffix  = defined $config->{output_suffix} ? $config->{output_suffix} : "filtered";
    my $out     = defined $config->{output}        ? "$config->{output}"      : "$out_dir/$prefix\_$suffix\.bam";
    my $cmd     = "samtools view -H $input";
    my $header  = $self->_run_cmd($cmd);
    open( my $in, "-|", "samtools view $input" )
        or $self->fail("&filter_bam_by_ids can't open input bam: $input because: $!\n");
    open( my $fh, "| samtools view -S - -bo $out" )
        or $self->fail("&filter_bam_by_ids can't open  output bam: $out because: $!\n");
    print $fh "$header";

    while (<$in>) {
        chomp;
        my @fields = split(/\t/);
        if ( $config->{good_list} && $good_ids->{ $fields[0] } ) { print $fh "$_\n"; $found_ids{ $fields[0] }++; }
        if ( $config->{bad_list}  && !$bad_ids->{ $fields[0] } ) { print $fh "$_\n"; $found_ids{ $fields[0] }++; }
    }
    close $in;
    close $fh;
    my $count = 0;
    foreach my $keys ( keys %found_ids ) { $count++; }
    if ( $self->{verbose} ) { print STDERR "======== &filter_bam_by_ids: Finished ========\n"; }
    return {
        count => $count,
        bam   => $out
    };
}

=head2 &_read_ids

 Title   : _read_ids
 Usage   : my $id_hash = $self->_read_ids({list => <ID_LIST_FILE>})
 Function: Create a hash of id's from a list file
 Returns : A hash of ids
 Args    : A hash containing potentially several config options:
        list     => File with ids 
=cut

sub _read_ids {
    my ( $self, $config ) = @_;
    if ( !$config->{list} ) { $self->fail("Must pass &_read_ids a list =>\n"); }
    if ( $self->empty_chk( { input => $config->{list} } ) == 1 ) {
        carp "Warning: &_read_ids input: $config->{list} is empty\n";
    }
    my $file = $config->{list};
    my %hash;
    open IN, "<$file" or $self->fail("&_read_ids unable to open: $file because: $!\n");
    while (<IN>) {
        chomp;
        $hash{$_} = 1;
    }
    close IN;
    return \%hash;
}

=head2 &empty_chk

 Title   : empty_chk
 Usage   : last if($self->empty_chk({input => <FILE>}) == 1);
 Function: Check to make sure a file is not empty
 Returns : 1 = Empty; 0 = input is not empty
 Args    : A hash containing potentially several config options:
        input     => Check if input is empty 
=cut

sub empty_chk {
    my ( $self, $config ) = @_;
    if ( !$config->{input} ) { $self->fail("Must pass &empty_chk a file."); }
    my $file  = $config->{input};
    my $empty = 0;                  ## 0 = False, 1=True, file is empty.
    my $count;
    if ( $file =~ /\.bam$/ ) {
        $count = `samtools view $file | head | wc -l`;
    }
    elsif ( $file =~ /\.gz/ ) {
        $count = `zcat $file | head | wc -l`;
    }
    else {
        $count = `head $file | wc -l`;
    }
    if ($?) { $self->fail("&empty_chk could not check: $config->{input}\n"); }
    if ( $count == 0 ) { $empty = 1; }
    return $empty;
}

=head2 &fail

 Title   : fail
 Usage   : if(1+1!=2){fail("Because you can't add!");}
 Function: Die with either die or confess depending on verbosity.
 Args    : Die message.
=cut

sub fail {
    my ( $self, $message ) = @_;
    chomp($message);
    if ( $self->{verbose} )      { confess "$message\n"; }
    if ( $self->{verbose} != 1 ) { die "$message\n"; }
}

=head2 &validated_bam

 Title   : validated_bam
 Usage   : my $validated_bam = lgtseek->validated_bam({by_clone=> <input> });
 Function: Creates a bam based on LGTFinder Results for reads with Human & Bacteria LGT
 Returns:  Hash->{count} && Hash->{file} bam based on LGTFinder output. 
 Args    : A hash containing potentially several config options:
        input           =>  Input bam to parse LGT reads from
        by_clone        =>  LGTFinder Output
        by_trace        =>  Not implemented yet.
        output_dir      =>  
        output_prefix   =>
        output_suffix   =>
        output          =>
=cut

sub validated_bam {
    my ( $self, $config ) = @_;
    if ( !$config->{input} ) {
        $self->fail("Must pass &validated_bam an input bam to parse reads from with input =>.\n");
    }
    if ( !$config->{by_clone} && !$config->{by_trace} ) {
        $self->fail("Must pass &validated_bam an LGTFinder output with by_clone => or by_trace=>.\n");
    }
    my $input = "$config->{input}";
    my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
    my $out_dir = defined $config->{output_dir}    ? $config->{output_dir}    : $path;
    my $prefix  = defined $config->{output_prefix} ? $config->{output_prefix} : $fn;
    my $suffix  = defined $config->{output_suffix} ? $config->{output_suffix} : "filtered";
    my $out     = defined $config->{output}        ? "$config->{output}"      : "$out_dir/$prefix\_$suffix\.bam";
    if ( $config->{output} ) {
        my ( $foo, $bar, $salad ) = fileparse( $config->{output}, qr/\.[^\.]+/ );
        $out_dir = $bar;
    }

    open( IN, "<", "$config->{by_clone}" )
        || $self->fail("ERROR: &validated_bam can not open input: $config->{by_clone}\n");
    open( OUT, ">", "$out_dir/$prefix\_valid_blast_read.list" )
        || $self->fail("ERROR: &validated_bam can not open output: $out_dir/$prefix\_valid_blast.list\n");
    ##  Make a list of reads from lgt_finder that have a "valid blast."
    while (<IN>) {
        chomp;
        my ( $read, $otu1, $otu2, $lca1, $lca2 ) = ( split /\t/, $_ )[ 0, 1, 2, 6, 11 ];

        # This needs to be checked. Might be too stringent? 12.04.13 KBS. This is also super specific for human lgt ...
        next if ( $otu1 =~ /Homo/      && $otu2 =~ /Homo/ );
        next if ( $otu1 !~ /Homo/      && $otu2 !~ /Homo/ );
        next if ( $lca1 =~ /Eukaryota/ && $lca2 =~ /Eukaryota/ );
        next if ( $lca1 =~ /Bacteria/  && $lca2 =~ /Bacteria/ );
        print OUT "$read\n";
    }

    close IN || $self->fail("ERROR: &validated_bam can't close input by_clone: $config->{by_clone} because: $!\n");
    close OUT
        || $self->fail(
        "ERROR: &validated_bam can't close output valid_blast.list: $out_dir/$prefix\_valid_blast_read.list because: $!\n"
        );

    ## Use the list of blast validated LGT's create a new bam.
    ## $valid_bam is a hash->{count} & ->{file}
    my $valid_bam = $self->filter_bam_by_ids(
        {   input_bam => $input,
            good_list => "$out_dir/$prefix\_valid_blast_read.list",
            output    => $out
        }
    );

    $self->_run_cmd("rm $out_dir/$prefix\_valid_blast_read.list");
    return $valid_bam;
}

=head2 new2

 Title   : new2
 Usage   : my $lgtseek = LGTSeek->new2(\%options)
 Function: Creates a new LGTSeek object, key=>values take #1= %options, #2=  ~/.lgtseek.config 
 Returns : An instance of LGTSeek

=cut

sub new2 {
    my ( $class, $options ) = @_;
    use File::HomeDir;

    # Usefull list for fileparse: @{$lgtseek->{'list'}}
    my $self = {
        bam_suffix_list => [
            '_resorted.\d+.bam', '_resorted.bam', '.gpg.bam',   '_prelim.bam',
            '_name-sort.bam',    '_pos-sort.bam', '_psort.bam', '_nsort.bam',
            '.srt.bam',          '.bam'
        ],
        sam_suffix_list   => [ '.sam.gz',       '.sam' ],
        fastq_suffix_list => [ '_\d+.fastq.gz', '.fastq.gz', '_\d+.fastq', '.fastq', '.fq' ],
        fasta_suffix_list   => [ '.fa.gz',   '.fasta.gz',     '.fasta', '.fa' ],
        mpileup_suffix_list => [ '.mpileup', '_COVERAGE.txt', '.txt' ],
        suffix_regex        => qr/\.[^\.]+/,
    };

    ## Now open the config file
    ## First determine the proper file path for the config file
    if ( $options->{conf_file} =~ /^\~(.*)/ ) {
        $options->{conf_file} = File::HomeDir->my_home . $1;
    }    ## This is incase the user passed a config file like ~/config.file
    my $conf_file
        = defined $options->{conf_file}
        ? $options->{conf_file}
        : File::HomeDir->my_home . "/.lgtseek.conf";    ## Open --conf_file or the default ~/.lgtseek.conf
    ## Open the config file and build a hash of key=>value for each line delimited on white space
    if ( -e $conf_file ) {
        my %config;
        open( IN, "<", "$conf_file" ) or confess "Can't open conf_file: $conf_file\n";
        while (<IN>) {
            chomp;
            next if ( $_ =~ /^#/ );
            my ( $key, $value ) = split( /\s+/, $_ );
            $config{$key} = $value;
        }
        close IN or confess "Error: can't close conf_file: $conf_file\n";
        ## Make sure all keys from --options and config.file have a value, priority goes to --option
        ## If a key was passed from --options use the --option=>value if avail, or use config.file=>value
        foreach my $opt_key ( keys %$options ) {
            $self->{$opt_key} = defined $options->{$opt_key} ? $options->{$opt_key} : $config{$opt_key};
        }
        ## Make sure all the keys from the config file have a value, use --options=>value if avail or use the config.file=>value
        foreach my $conf_key ( keys %config ) {
            $self->{$conf_key} = defined $options->{$conf_key} ? $options->{$conf_key} : $config{$conf_key};
        }
    }
    else {
        $self = $options;
    }

    bless $self;
    return $self;
}

1;

