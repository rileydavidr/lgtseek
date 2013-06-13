#!/usr/bin/perl

=head1 NAME

start_lgtmongo.pl

=head1 SYNOPSIS

Starts a mongo instance on a clovr host and loads a gi2taxon database

=head1 DESCRIPTION

Use this script to start/load a gi2taxon database. It assumes you have started 
a cluster.

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
                          'taxon_host=s', # Comma separated list of files
                          'taxon_dir=s',
                          'port=s',
                          'taxon_idx_dir=s'
                          );

my $port=$options{port} ? $options{port} : 10000;

# First prevent this machine from automatically shutting down
print `ssh -oNoneSwitch=yes -oNoneEnabled=yes -o PasswordAuthentication=no -o ConnectTimeout=30 -o StrictHostKeyChecking=no -o ServerAliveInterval=30 -o UserKnownHostsFile=/dev/null -q -i /mnt/keys/devel1.pem root\@$options{taxon_host} 'touch /var/vappio/runtime/noautoshutdown'`;

# First start a mongo instance on the remote host.
print `ssh -oNoneSwitch=yes -oNoneEnabled=yes -o PasswordAuthentication=no -o ConnectTimeout=30 -o StrictHostKeyChecking=no -o ServerAliveInterval=30 -o UserKnownHostsFile=/dev/null -q -i /mnt/keys/devel1.pem root\@$options{taxon_host} \'mkdir /mnt/lgtmongo;mongod --quiet --dbpath=/mnt/lgtmongo --logpath=/mnt/lgtmongo.log --fork --port=$port\'`;
if($?) {
    die "$?";
}
`sleep 10`;
`mongo $options{taxon_host} -eval "db.currentOp()" --port $port`;

# We'll try again if we were unable on the last one. If this doesn't work we give up.
if($?) {
    print `ssh -oNoneSwitch=yes -oNoneEnabled=yes -o PasswordAuthentication=no -o ConnectTimeout=30 -o StrictHostKeyChecking=no -o ServerAliveInterval=30 -o UserKnownHostsFile=/dev/null -q -i /mnt/keys/devel1.pem root\@$options{taxon_host} \'mkdir /mnt/lgtmongo;mongod --quiet --dbpath=/mnt/lgtmongo --logpath=/mnt/lgtmongo.log --fork --port=$port\'`
}

# Create an lgtseek object
my $lgtseek = LGTSeek->new({
    taxon_dir => $options{taxon_dir},
    taxon_host => $options{taxon_host}.":$port",
    chunk_size => '500000',
    taxon_idx_dir => $options{taxon_idx_dir}, 
});

my $gi2tax = $lgtseek->getGiTaxon();
# OK, hopefully that worked.
