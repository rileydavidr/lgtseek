package setup_default_paths;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( setup_default_paths );
## This module is handy for setting default paths for scripts using lgt_seek.pm
## Assumes the main script uses $options{fs || diag || clovr} = <0|1> [0]

sub setup_default_paths {
	our %options;
	our $results;
	my $diag = 
	{
		bin_dir => "/opt/lgtseek/bin/",
		ergatis_bin => "/opt/ergatis/bin/",
		prinseq_bin => "/opt/prinseq/bin/",
		samtools_bin => "samtools",
		split_bac_list => "/mnt/staging/data/lgt_seq/mnt/references/split_refseq_bacteria/split_bacteria_ref.list",
		hg19_ref => "/mnt/staging/data/lgt_seq/mnt/references/hg19/dna/hg19.fa",
		refseq_list => "/mnt/staging/data/lgt_seq/mnt/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "cloud-128-152.diagcomputing.org:10001",
		taxon_dir => "/mnt/staging/data/lgt_seq/mnt/references/taxonomy",
		taxon_idx_dir => "/mnt/staging/data/lgt_seq/mnt/references/taxonomy",
		path_to_blastdb => "/mnt/staging/data/ncbi-nt/nt",
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};
	
	my $clovr = 
	{
		bin_dir => "/opt/lgtseek/bin/",
		ergatis_bin => "/opt/ergatis/bin/",
		prinseq_bin => "/opt/prinseq/bin/",
		samtools_bin => "samtools",
		split_bac_list => "/mnt/references/split_refseq_bacteria/split_bacteria_ref.list",
		hg19_ref => "/mnt/references/hg19/dna/hg19.fa",
		refseq_list => "/mnt/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "cloud-128-152.diagcomputing.org:10001",						
		taxon_dir => "/mnt/references/taxonomy/20120720_135302/",
		taxon_idx_dir => "/mnt/references/taxonomy/20120720_135302/",
		path_to_blastdb => "/mnt/scratch/ksieber/ref/refseq_bacteria_merged.fna",  ## FIX this
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};

	my $fs = 
	{
		bin_dir => "/local/projects-t3/HLGT/scripts/lgtseek/bin/",
		ergatis_bin => "/local/projects/ergatis/package-driley/bin/",
		prinseq_bin => "/home/ksieber/lib/prinseq-lite-0.18.1/",
		samtools_bin => "samtools",
		split_bac_list => "/local/projects-t3/HLGT/references/split_bacteria/all_bacteria.list",
		hg19_ref => "/local/projects-t3/HLGT/references/hg19/hg19.fa",
		refseq_list => "/local/projects-t3/HLGT/references/refseq_bacteria_BWA_INDEXED_20110831/refseq.list",
		taxon_host => "mongotest1-lx.igs.umaryland.edu:10001",
		taxon_dir => "/local/db/repository/ncbi/blast/20120414_001321/taxonomy/",
		taxon_idx_dir => "/local/projects-t3/HLGT/idx_dir/20120414",
		path_to_blastdb => "/local/db/ncbi/blast/db/nt",  ## FIX this
		donor_lineage => "Bacteria",
		host_lineage => "Eukaryota"
	};
	
	if($main::options{clovr}==1){foreach my $keys (keys %$clovr){$main::options{$keys}=$clovr->{$keys};}}
	if($main::options{diag}==1){foreach my $keys (keys %$diag){$main::options{$keys}=$diag->{$keys};}}			
	if($main::options{fs}==1){foreach my $keys (keys %$fs){$main::options{$keys}=$fs->{$keys};}}
}
1;
