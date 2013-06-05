#!/usr/bin/perl

=head1 NAME

lgtseek_microbiome_lca.pl

=head1 SYNOPSIS

Takes a g-zipped sam file and runs BWA against specified references. Filters with PRINSeq and reads mapping to another reference, then calculates the LCA.

=head1 DESCRIPTION

Takes a g-zipped sam file. Runs BWA searches against bacteria in RefSeq and the human genome and then removes read pairs where either read align to the human genome. Also filters low complexity and duplicate reads with PRINSeq, then calculates the LCA.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use strict;
use LGTSeek;
use lib '../lib';
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
#                          'run_id=s',
                          'input_files=s',   # Comma separated list of sam.gz files
                          'donor_files=s',   # Directory for references in RefSeq
                          'host_file=s',     # Path to hg19 file for BWA alignment
                          'output_dir=s',    # Output directory
                          'bin_dir=s',
                          'sratools_dir=s',
                          'samtools_bin=s',
                          'aspera_dir=s',
                          'ergatis_dir=s',
                          'prinseq_bin=s',
                          'help|h'
                          );


# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects/ergatis/package-driley/bin/'; 

my $ergatis_dir = $options{ergatis_dir} ? $options{ergatis_dir} :'/local/projects/ergatis/package-driley/bin/';

my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : 'prinseq-lite.pl';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    ergatis_bin => $ergatis_dir,
    prinseq_bin => $prinseq_bin,
    paired_end => 1
});

#my @refs = split(',',$options{ref_files}); #I was going to change the names of the arrays to make it easier to understand what I'm doing, but I think they need to have the same names in order for the subroutines to run correctly
#my @filter_refs = split(',',$options{filter_ref_files});

my @donor_refs = split(',',$options{donor_files});
my @host_refs = split(',',$options{host_files});

# Now we'll loop over the run files (in case there are multiple)
foreach my $file (@$input_files) {

    my($name,$path,$suff) = fileparse($file,'.sam.gz');
    chomp $name;
    # Next convert from '.sam.gz' to .bam
    my $bam = "$path/$name.bam";
    `samtools view -bS -o $bam $file`;
    ########################## Need to come back and finish this part, it would currently work if the sam.gz files had headers, but none of the ones that I have in the KIPR_seq directory have them, so I can't make them into bam files. Is there a quick way to get around this or to add a header to all of the files? Or will the files on diag be differet? #################################

# It might make more sense to prinseq and filter for human before running the runBWA with all of RefSeq, this way we won't have to filter out each individual bam file from RefSeq alignments. It shouldn't be too hard to do the prinseq from bam right after the file is converted to a bam file, and then map to refseq and human like I have it, then choose the MM_UU files

    # Filter with PRINSeq for low complexity and duplicates.
    my $prinseq_bams = $lgtseek->prinseqFilterBam(
        {bam_file => $bam,
         output_dir => "$options{output_dir}/prinseq_filt_bams/"
        });

    print STDERR "Prinseq filtered bams: @$prinseq_bams\n";

    # Align to the human reference.
    my $host_bams = $lgtseek->runBWA(
        {input_bam => "$options{output_dir}/prinseq_filt_bams/$name\_filtered.bam", # this should be the output from the prinseq filtering
         output_dir => "$options{output_dir}/host_alignments/",
         output_bam => 1,
         other_opts => {threads => 6},
         reference => join(',',@host_refs),
         overwrite => 1
        });

    print STDERR "Host bams: @$host_bams\n";

   # Filter for reads that mapped to human.
    my $filtered_bam = $lgtseek->filterBamfromBam(
        {output_dir=> "$options{output_dir}/filtered_bam/",
         filter => "$options{filter}",
         input => "$options{input}"
        });

    print STDERR "Filtered bams: @$filtered_bams\n";

    # Align to the RefSeq.
    my $donor_bams = $lgtseek->runBWA(
        {input_bam => $bam, 
         output_bam => 1,
         other_opts => {threads => 4},
         output_dir => "$options{output_dir}/donor_alignments/",
         reference => join(',',@donor_refs)
        });

    print STDERR "RefSeq bams: @$bams\n";

    # Postprocess the results
    my $pp_data = $lgtseek->bwaPostProcess(
        {donor_bams => $donor_bams,
         host_bams => $host_bams,
         output_prefix => $name,
         overwrite => 1
        });

    print STDERR "Postprocess: @$pp_data\n";

