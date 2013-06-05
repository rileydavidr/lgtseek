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
use lib '../lib';
use strict;
use LGTSeek;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions (\%options,
                          'bam_file=s', # Comma separated list of files
                          'ref_db=s',
                          'output_dir=s',
                          'bin_dir=s',
                          'samtools_bin=s',
                          'output_dir=s',
                          'donor_lineage=s',
                          'host_lineage=s',
                          'taxon_host=s',
                          'taxon_dir=s',
                          'taxon_idx_dir=s',
                          'help|h'
                          );

# Take care of the inputs
my $bin_dir = $options{bin_dir} ? $options{bin_dir} : '/local/projects/ergatis/package-driley/bin/';

my $sra_dir = $options{sra_dir} ? $options{sra_dir} : '/usr/local/packages/sratoolkit.2.1.8/';

my $aspera_dir = $options{aspera_dir} ? $options{aspera_dir} : '/home/driley/.aspera/connect/';
    
my $ergatis_dir = $options{ergatis_dir} ? $options{ergatis_dir} :'/local/projects/ergatis/package-driley/bin/';

my $prinseq_bin = $options{prinseq_bin} ? $options{prinseq_bin} : 'prinseq-lite.pl';

my $samtools_bin = $options{samtools_bin} ? $options{samtools_bin} : 'samtools';




# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    bin_dir => $bin_dir,
    output_dir => $options{output_dir},
    aspera_path => $aspera_dir,
    ergatis_bin => $ergatis_dir,
    prinseq_bin => $prinseq_bin,
    paired_end => 1,
    samtools_bin => $samtools_bin,
    taxon_host => $options{taxon_host},
    taxon_dir => $options{taxon_dir},
    taxon_idx_dir => $options{taxon_idx_dir}
});

my($name,$path,$suff) = fileparse($options{bam_file},'.bam');
# Process the LGT's here
my $filtered_fasta = $lgtseek->sam2Fasta({
    input => $options{bam_file}});
    
# Now blast validate the putative LGTs

# First get the best hits
my $best_blasts = $lgtseek->bestBlast2(
    {
        db => $options{path_to_nt},
        lineage1 => $options{donor_lineage},
        lineage2 => $options{host_lineage},
        fasta => $filtered_fasta,
        output_dir => "$options{output_dir}"
    });

# Now run lgtfinder
my $valid_lgts = $lgtseek->runLgtFinder(
    {
        lineage1 => 'Eukaryota',
        lineage2 => 'Bacteria',
        input_file_list => $best_blasts->{list_file},
        output_prefix => "$name"
    });

# Also run the blast again with raw output.
`blastall -p blastn -e 10e-5 -T F -d $options{path_to_nt} -i $filtered_fasta > $options{output_dir}/$name\_blast.raw`;



#push(@header,'valid_lgt');
#push(@vals,$valid_lgts->{valid_clones});

#&print_tab("$options{output_dir}/$name\_blast_validate.tab",\@header,\@vals);

#sub print_tab {
#    my ($file,$header,$vals) = @_;
#    open OUT, ">$file" or die "Couldn't open $file\n";
#    print OUT join("\t",@$header);
#    print OUT "\n";
#    print OUT join("\t",@$vals);
#    print OUT "\n";    
#}
