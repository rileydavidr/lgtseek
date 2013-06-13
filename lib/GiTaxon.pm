=head1 NAME

GiTaxon

=head1 DESCRIPTION

Utility to look up taxon information in a mongo database.

=head1 AUTHOR

David Riley

driley@som.umaryland.edu

=cut

package GiTaxon;
use strict;
use MongoDB;
use Bio::DB::Taxonomy;
use Bio::DB::EUtilities;
use File::Find;

sub new {
    my ($class, $args) = @_;

    my $self = {};
    $self->{'nodes'} = $args->{'nodes'} ? $args->{'nodes'} : '/local/db/by_source/ncbi/taxonomy/latest/names.dmp';
    $self->{'names'} = $args->{'names'} ? $args->{'names'} : '/local/db/by_source/ncbi/taxonomy/latest/names.dmp';
    $self->{'gi2tax'} = $args->{'gi2tax'} ? $args->{'gi2tax'} : '/local/db/by_source/ncbi/taxonomy/latest/gi_taxid_nucl.dmp';
    $self->{'chunk_size'} = $args->{'chunk_size'} ? $args->{'chunk_size'} : 10000;
    $self->{'idx_dir'} = $args->{'idx_dir'} ? $args->{'idx_dir'} : '/tmp/';
    $self->{'host'} = $args->{'host'} ? $args->{'host'} : 'tettelin-lx.igs.umaryland.edu';
    $self->{'gi_db'} = $args->{'gi_db'} ? $args->{'gi_db'} : 'gi2taxon';
    $self->{'gi_coll'} = $args->{'gi_coll'} ? $args->{'gi_coll'} : 'gi2taxonnuc';
    $self->{'taxonomy_dir'} = $args->{'idx_dir'} ? $args->{'idx_dir'} : '/tmp';

    $self->{'cache'} = {}; 

    $self->{'type'} = $args->{'type'} ? $args->{'type'} : 'nucleotide';
    my $gi_tax_file = 'gi_taxid_nucl.dmp';
    if($self->{'type'} eq 'protein') {
        $gi_tax_file = 'gi_taxid_prot.dmp';
    }
    print STDERR "Here with $args->{'taxon_dir'}\n";
    # This option can be used if the user want's to override all the nodes/names params at once
    if($args->{'taxon_dir'}) {
        print STDERR "Here with a taxon directory $args->{'taxon_dir'}\n";
        # Find the nodes, names and nucleotide mapping file
        find(sub {
            if($File::Find::name =~  /nodes.dmp/) {
                
                $self->{'nodes'} = $File::Find::name;
            }
            elsif($File::Find::name =~ /names.dmp/) {
                $self->{'names'} = $File::Find::name;
            }
            elsif($File::Find::name =~ /$gi_tax_file/) {
                $self->{'gi2tax'} = $File::Find::name;
            }},$args->{'taxon_dir'});
        if(!$args->{'gi_coll'}) {
            if($args->{'taxon_dir'} =~ /(\d+\_\d+)/) {
                my $date = $1;
                if($self->{'type'} eq 'nucleotide') {
                    $self->{'gi_coll'} = "gi2taxonnuc_$date";
                }    
                else {
                    $self->{'gi_coll'} = "gi2taxonprot_$date";
                }
            }    
        }
    }
    if($args->{'nodes'}) {
        $self->{'nodes'} = $args->{'nodes'};
    }
    if($args->{'names'}) {
        $self->{'names'} = $args->{'names'};
    }
    if($args->{'gi2tax'}) {
        $self->{'gi2tax'} = $args->{'gi2tax'};
    }
    $self->{'db'} = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                 -nodesfile => $self->{'nodes'},
                                 -namesfile => $self->{'names'},
                                 -directory => $self->{'taxonomy_dir'});
    print STDERR "Using $self->{nodes} and $self->{names} and $self->{'taxonomy_dir'} $self->{'gi_coll'} $self->{'db_host'} $self->{'gi2tax'}\n";
    bless $self;
    $self->{'gi2taxon'} = $self->getgi2taxon($self->{'gi2tax'});

    return $self;
}

sub getTaxon {
    my($self,$acc) = @_;
    $acc =~ s/^\s+//;
    $acc =~ s/\s/\t/;
    my $gi = $acc;
    if($acc =~ /\|/) {
        my @fields = split(/\|/, $acc);
        $gi = $fields[1]; 
    }
    my $taxonid = '';
    my $retval = {};

    # First check the cache
    if($self->{cache}->{$gi}) {
        $retval = $self->{cache}->{$gi};
    }
    else {
        my $taxon_lookup = $self->{'gi2taxon'}->find_one({'gi' => "$gi"}, {'taxon' => 1});
        
        
        if($taxon_lookup) {
            $taxonid = $taxon_lookup->{'taxon'};
        }
        else {
            print STDERR "Unable to find taxon for $gi, Checking NCBI\n";
            my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
#                                       -email => 'mymail@foo.bar',
                                       -db    => $self->{'type'},
                                       -id    => [$gi]);
            while (my $ds = $factory->next_DocSum) {
                my @res = $ds->get_contents_by_name('TaxId');
                if(@res) {
                    $taxonid = $res[0];
                }
                if(!$taxonid) {
                    print STDERR "Unable to find taxonid at NCBI\n";
                }
                else {
                    $self->{'gi2taxon'}->update({'gi' => "$gi"}, {'gi' => "$gi", 'taxon' => $taxonid}, {'upsert' => 1});
                    print STDERR "Added $gi\t$taxonid to the db\n";
                }
#                print "ID: ",$ds->get_id,"\n";
                # flattened mode
#                while (my $item = $ds->next_Item('flattened'))  {
#                    # not all Items have content, so need to check...
#                    printf("%-20s:%s\n",$item->get_name,$item->get_content) if $item->get_content;
#                }
#                print "\n";
            }

        }
        
        my $taxon = $self->{'db'}->get_taxon(-taxonid => $taxonid);
        if(!$taxon) {
#        print STDERR "Unable to find taxon for $gi\n";
        }
        if(!$taxon) {
            $retval = {
                'acc'=> $acc,
                'gi'=> $gi,
                'taxon_id'=> $taxonid
            };
        }
        elsif ($taxon->isa('Bio::Taxon')) {
            my $name = $taxon->scientific_name;
            my $c = $taxon;
            my @lineage = ($name);
            while (my $parent = $self->{'db'}->ancestor($c)) {
                unshift @lineage, $parent->scientific_name;
                $c = $parent;
            }
            $retval = {
                'gi' => $gi,
                'acc'=> $acc,
                'taxon_id'=> $taxonid, 
                'name'=> $name, 
                'lineage'=> join(";", @lineage)
            };
        }else {
            print STDERR "Had something other than a Bio::Taxon\n";
#	print join("\t", $taxonid),"\n";
        }
        $self->{cache}->{$gi} = $retval;
    }
    return $retval;
}
sub getgi2taxon {
    my($self, $data_file) = @_;

    my $mongo = $self->get_mongodb_connection($self->{'gi_db'},$self->{'host'});
    my $coll = $mongo->get_collection($self->{'gi_coll'});
    if(!$coll->find_one()) {
        print "Found nothing in database $self->{gi_db} collection $self->{gi_coll} on $self->{host}\n";
        print "Getting the line count\n";
        my $lc = `wc -l $data_file`;
        chomp $lc;
        $lc =~ s/\s.*//;
        print "Got the line count\n";
        open IN, "<$data_file" or die "Unable to open $data_file\n";
        my $num_in_chunk = 0;
        my $total = 0;
        my @chunk;
        while(<IN>) {
            chomp;
            my ($gi, $taxon) = split(/\t/, $_);
            $num_in_chunk++;
            push(@chunk, {'gi' => $gi, 'taxon' => $taxon});
            if($num_in_chunk==$self->{'chunk_size'}) {
                $total += $num_in_chunk;
                print join("", ("\r",(sprintf('%.2f', (($total/$lc)*100))),"% complete"));
                $self->insert_chunk($coll, \@chunk);
                @chunk = ();
                $num_in_chunk =0;
            }
        }
        $self->insert_chunk($coll,\@chunk);

        close IN;
        $coll->ensure_index({'gi' => 1},{'safe' => 1});
    }
    return $coll;
}
sub insert_chunk {
    my $self =shift;
    my $coll = shift;
    my $chunk = shift;
    $coll->batch_insert($chunk,{'safe' => 1});
}
sub get_mongodb_connection {
    my($self,$dbname,$host) = @_;
    # First we'll establish our connection to mongodb
    my $conn = MongoDB::Connection->new(host => $host);
    $conn->query_timeout(60000000);
    return $conn->get_database($dbname);
}

1;
