#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use JSON::XS;

my $filename = $ARGV[0];
open (my $fh,"<",$filename);
my $data;
while (my $line = <$fh>) {
	$data .= $line;	
}
close($fh);
my $JSON = JSON::XS->new->utf8(1);
$data = $JSON->decode($data);
