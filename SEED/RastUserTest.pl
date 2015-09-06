#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use File::Path;
use SOAP::Lite;
use Data::Dumper;

my $username = $ARGV[0];
my $dsn = "DBI:mysql:WebAppBackend:bio-app-authdb.mcs.anl.gov:3306";
my $user = "webappuser";
my $db = DBI->connect($dsn, $user);
if (!defined($db)) {
	die("Could not connect to database!");
}
my $select = "SELECT * FROM User WHERE User.login = ?";
my $columns = {};
my $users = $db->selectall_arrayref($select, { Slice => $columns }, $username);
$db->disconnect;
print Data::Dumper->Dump($users);