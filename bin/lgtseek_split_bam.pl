#!/usr/bin/perl

=head1 NAME

lgtseek_split_bam.pl

=head1 SYNOPSIS

Splits input bams into chunks of a certain size.

=head1 DESCRIPTION

Splits input bams into chunks of a certain size.

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
                          'bam_file_list=s',
                          'bam_files=s',
                          'output_dir=s',
                          'samtools_bin=s',
                          'ergatis_dir=s',
                          'output_list=s',
                          'help|h'
                          );

# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects/ergatis/package-driley/bin/';

my $sra_dir = $options{sra_dir} ? $options{sra_dir} : '/usr/local/packages/sratoolkit.2.1.8/';

my $ergatis_dir = $options{ergatis_dir} ? $options{ergatis_dir} :'/local/projects/ergatis/package-driley/bin/';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    ergatis_bin => $ergatis_dir,
    samtools_bin => $samtools_bin,
    paired_end => 1
});

# Deal with the input bams
my @bam_files = split(',',$options{bam_files});

if($options{bam_file_list}) {
    open IN, "<$options{bam_file_list}" or die "Couldn't open $options{bam_file_lsit}\n";
    while(<IN>) {
        chomp;
        push(@bam_files,$_);
    }
}

# Open a list file to write the output bams to
my $olistfh;
if($options{output_list}) {
    open($olistfh, "<$options{output_list}") or die "Unable to open $options{output_list}\n";
}

# Loop over the bam files and split them
foreach my $file (@bam_files) {
    my $split_files = $lgtseek->splitBam({
        input => $file,
        output_dir => $options{output_dir}
    });

    # Write out the filenames to the list file.
    if($olistfh) {
        map {
            print $olistfh "$_\n";
        }@$split_files
    }
}
