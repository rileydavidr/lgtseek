#!/usr/bin/perl

=head1 NAME

lgtseek_sam2lca.pl

=head1 SYNOPSIS

Map a bam file against a set of references and write out the lca's of all the
reads.

=head1 DESCRIPTION

Map a bam file against a set of references and write out the lca's of all the
reads.

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
use lib '../lib';
use strict;
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'input_bam=s', # Comma separated list of files
                          'reference_list=s',
                          'references=s',
                          'output_file=s',
                          'output_dir=s',
                          'samtools_bin=s',
                          'threads=s',
                          'taxon_host=s',
                          'taxon_dir=s',
                          'taxon_idx_dir=s',
                          'help|h'
                          );


# Take care of the inputs
my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : 'prinseq-lite.pl';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';

my $threads = $options{threads} ? $options{threads} : 1;

my @references;

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    output_dir => $options{output_dir},
    prinseq_bin => $prinseq_bin,
    taxon_host => $options{taxon_host},
    taxon_dir => $options{taxon_dir},
    taxon_idx_dir => $options{taxon_ids_dir},
    paired_end => 1
});


$lgtseek->runBWA({
    input_bam => $options{input_bam},
    out_file => $options{output_file},
    reference_list => $options{reference_list},
    cleanup_sai => 1,
    run_lca => 1
});

