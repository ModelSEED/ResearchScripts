#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];

my $sapsvr = ModelSEED::Client::SAP->new();
my $genomeHash = $sapsvr->all_genomes({});
my $list = [keys(%{$genomeHash})];
open (my $out, ">", $directory."/PubSEED.txt");
for (my $i=0; $i < @{$list}; $i++) {
	print $out $list->[$i]."\t".$genomeHash->{$list->[$i]}."\n";
}
close($out);