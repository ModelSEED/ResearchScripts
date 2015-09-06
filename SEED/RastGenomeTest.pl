#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use File::Path;
use SOAP::Lite;
use Data::Dumper;

my $genome = $ARGV[0];
my $dsn = "DBI:mysql:RastProdJobCache:rast.mcs.anl.gov:3306";
my $user = "rast";
my $db = DBI->connect($dsn, $user);
if (!defined($db)) {
	die("Could not connect to database!");
}
my $select = "SELECT * FROM Job WHERE Job.genome_id = ?";
my $columns = {
#	_id		 => 1,
#	id		  => 1,
#	genome_id   => 1
};
my $jobs = $db->selectall_arrayref($select, { Slice => $columns }, $genome);
$db->disconnect;
print Data::Dumper->Dump($jobs);