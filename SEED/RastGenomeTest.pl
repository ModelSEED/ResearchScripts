#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use File::Path;
use SOAP::Lite;
use Data::Dumper;

my $username = $ARGV[0];
my $password = $ARGV[1];

my $seed = join '', ('.', '/', 0..9, 'A'..'Z', 'a'..'z')[rand 64, rand 64];
print crypt($password, $seed);

my $dsn = "DBI:mysql:WebAppBackend:bio-app-authdb.mcs.anl.gov:3306";
my $user = "webappuser";
my $db = DBI->connect($dsn, $user);
if (!defined($db)) {
	die("Could not connect to database!");
}
my $select = "SELECT * FROM User WHERE User.login = ?";
my $columns = {
	_id       => 1,
	login     => 1,
	password  => 1,
	firstname => 1,
	lastname  => 1,
	email     => 1
};
my $users = $db->selectall_arrayref($select, { Slice => $columns }, $username);
if (!defined($users) || scalar @$users == 0) {
	die("Username not found!");
}
$db->disconnect;
print Data::Dumper->Dump($users);