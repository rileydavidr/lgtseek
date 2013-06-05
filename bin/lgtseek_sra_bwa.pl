#!/usr/bin/perl

=head1 NAME

lgtseek_sra_bwa.pl

=head1 SYNOPSIS

Downloads an SRA file from the Sequence read archive and runs a bwa-based
LGT finding method.

=head1 DESCRIPTION

Downloads the SRA file(s) associated with the run_id provided. Runs
BWA searches against a set of donor and host references and then categorizes
each read as lgt,integration site, microbiome, host, no_map, all_map and single_map.

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use strict;
use lib '../lib';
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'run_id=s', # Comma separated list of files
                          'donor_files=s',
                          'host_files=s',
                          'output_dir=s',
                          'bin_dir=s',
                          'sratools_dir=s',
                          'samtools_bin=s',
                          'aspera_dir=s',
                          'ergatis_dir=s',
                          'prinseq_bin=s',
                          'threads=s',
                          'help|h'
                          );


# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects/ergatis/package-driley/bin/';

my $sra_dir = $options{sratools_dir} ? $options{sratools_dir} : '/usr/local/packages/sratoolkit.2.1.8/';

my $aspera_dir = $options{aspera_dir} ? $options{aspera_dir} : '/home/driley/.aspera/connect/';
    
my $ergatis_dir = $options{ergatis_dir} ? $options{ergatis_dir} :'/local/projects/ergatis/package-driley/bin/';

my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : 'prinseq-lite.pl';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';

my $threads = $options{threads} ? $options{threads} : 1;

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    sratoolkit_path => $sra_dir,
    aspera_path => $aspera_dir,
    ergatis_bin => $ergatis_dir,
    prinseq_bin => $prinseq_bin,
    paired_end => 1
});

my @donor_refs = split(',',$options{donor_files});
my @host_refs = split(',',$options{host_files});

# First download SRA files
my $sra_files = $lgtseek->downloadSRA(
    {run_id => $options{run_id}});

# Now we'll loop over the run files (in case there are multiple)
foreach my $file (@$sra_files) {

    my($name,$path,$suff) = fileparse($file,'.sra');
    chomp $name;
    # Next dump fastqs
    my $fastq_info = $lgtseek->dumpFastq(
        {sra_file => $file});

    if(!$fastq_info->{paired_end}) {
        print STDERR "Had an unpaired file $file\n";
        next;
    }
    

    # Align to the donors.
    my $donor_bams = $lgtseek->runBWA(
        {input_dir => $fastq_info->{path},
         input_base => $fastq_info->{basename},
         output_bam => 1,
         threads => $threads,
         output_dir => "$options{output_dir}/donor_alignments/",
         reference => join(',',@donor_refs)
        });

    print STDERR "Donor bams: @$donor_bams\n";
    # Align to the hosts.
    my $host_bams = $lgtseek->runBWA(
        {input_dir => $fastq_info->{path},
         input_base => $fastq_info->{basename},
         threads => $threads,
         output_dir => "$options{output_dir}/host_alignments/",
         output_bam => 1,
         reference => join(',',@host_refs)
        });
    print STDERR "Host bams: @$host_bams\n";

    # Postprocess the results
    my $pp_data = $lgtseek->bwaPostProcess(
        {donor_bams => $donor_bams,
         host_bams => $host_bams,
         output_prefix => $name,
         overwrite => 1
        });

    print STDERR "Removing the raw donor/host mappings\n";
    print STDERR `rm -rf $options{output_dir}/host_alignments/`;
    print STDERR `rm -rf $options{output_dir}/host_alignments/`;

    my @header = ('run_id');
    my @vals = ($name);
    open OUT, ">$options{output_dir}/$name\_post_processing.tab" or die;
    map {
        push(@header,$_);
        push(@vals,$pp_data->{counts}->{$_});
    } ('total','host','no_map','all_map','single_map','integration_site_host','integration_site_donor','microbiome','lgt');
    &print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

    # Prinseq filter the putative lgts
    my $filtered_bam = $lgtseek->prinseqFilterBam(
        {output_dir => "$options{output_dir}/prinseq_filtering",
         bam_file => $pp_data->{files}->{lgt_donor}});
    # Add filtered count to counts.
    push(@header,'prinseq_filtered_lgt');
    push(@vals,$filtered_bam->{count});
    &print_tab("$options{output_dir}/$name\_post_processing.tab",\@header,\@vals);

}

sub print_tab {
    my ($file,$header,$vals) = @_;
    open OUT, ">$file" or die "Couldn't open $file\n";
    print OUT join("\t",@$header);
    print OUT "\n";
    print OUT join("\t",@$vals);
    print OUT "\n";
}

