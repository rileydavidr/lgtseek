#!/usr/bin/perl

=head1 NAME

lgtseek_downloadCGHub.pl

=head1 SYNOPSIS

Takes a list of ID's and downloads the corresponding bam and bai files from TCGA's CGHub.

=head1 DESCRIPTION

Takes a list of either analysis ID's, URI's, or gto files, or an xml file from cgquery with everything you want to download. It also takes a path to the cghub.key file that provides the user's credentials to CGHub. The user must also provide it with an output directory and in some cases a path to the GeneTorrent installation.  

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use strict;
use lib '/opt/lgtseek/lib';
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'download_list=s',   # List of ID's with each ID on a newline.
                          'output_dir=s',      # Output directory
                          'bin_dir=s',
                          'genetorrent_path=s',
                          'cghub_key=s',
                          'help|h'
                          );


# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/opt/lgtseek/bin/';

my $gtpath = $options{genetorrent_path} ? $options{genetorrent_path} : '/opt/genetorrent/bin';

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    paired_end => 1,
    genetorrent_path => $gtpath
});

my $downloads = $options{download_list};
open (IN, "<", $downloads) or 
    die "Unable to open $downloads for reading.\n";

# Now we'll loop over the download ID's (in case there are multiple).
while (my $id = <IN>) {
    chomp $id;
    print "$id\n";
    # Download the files for that ID from CGHub.
    my $files = $lgtseek->downloadCGHub(
        {cghub_key => $options{cghub_key},
         analysis_id => $id,
         output_dir => $options{output_dir}
        });
    print STDERR "Downloaded bam files: ".@{$files->{bam_files}}."\n";
    print STDERR "Downloaded bai files: ".@{$files->{bai_files}}."\n";
    print STDERR "Downloaded gto files: ".@{$files->{gto_files}}."\n";
}
