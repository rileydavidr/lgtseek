#!/usr/bin/perl


=head1  NAME 

filter_sam_from_prinseq.pl

=head1 SYNOPSIS

USAGE: filter_sam_from_prinseq.pl 
        --good_list=path_to_good_list
        --bad_list=path_to_bad_list
        --sam_list=path_to_sam_list
        --sam_file=path_to_sam_file
        --out_file=path_to_out_file

=head1 OPTIONS

B<--accession_list> 
    Path to an accession list

B<--help,-h> 
    This help message/documentation.

=head1   DESCRIPTION

not written

=head1 CONTACT

    David Riley
    daveriley@users.sf.net

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
$|++;

my %options = ();
my $results = GetOptions (\%options,
              'good_list:s',
              'bad_list:s',
              'sam_list:s',
              'sam_file:s',
              'out_file:s',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $good_ids = {};
my $bad_ids = {};
my $samtools = 'samtools';
my $out = \*STDOUT;
if($options{out_file}) {
    open($out, ">$options{out_file}") or die "Unable to open $options{out_file}\n";
}

if($options{good_list}) {
    &read_ids($good_ids,$options{good_list});
}

if($options{bad_list}) {
    &read_ids($bad_ids,$options{bad_list});
}

my @files;
if($options{sam_list}) {
    open LIST, "<$options{sam_list}" or die "Unable to open list $options{sam_list}\n";
    @files = <LIST>;
    close LIST;
}
if($options{sam_file}) {
    push(@files, $options{sam_file});
}

foreach my $file (@files) {
    chomp $file;
    &read_sam($file);
}

sub read_sam {
    my $file = shift;

    my $handle;
    if($file =~ /.bam/) {
        open($handle, '-|', "$samtools view -h $file") or die "Unable to open $file with $samtools\n";
    }
    else {
        open($handle, "<$file") or die "Unable to open $file\n";
    }
    while(<$handle>) {
        chomp;
        my @lines = split(/\t/);
        if(/^@/) {
            print "$_\n";
        }
        if($options{good_list} && $good_ids->{$lines[0]}) {
            print "$_\n";
        }
        elsif($options{bad_list} && !$bad_ids->{$lines[0]}) {
            print "$_\n";
        }
    }
}

sub read_ids {
    my( $hash, $file ) = @_;

    open IN, "<$file" or die "Unable to open $file\n";

    while(<IN>) {
        chomp;
        $hash->{$_} = 1;
    }
    close IN;
}


