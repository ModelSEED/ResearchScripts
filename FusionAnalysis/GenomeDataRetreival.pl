#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];

my $sapsvr = ModelSEED::Client::SAP->new();
my $data = my $genomeHash = $sapsvr->all_genomes({});
my $list = [keys(%{$data})];
my $genomeHash = $sapsvr->genome_data({
	-ids => $list,
	-data => ["pegs"]
});
open (my $out, ">", $directory."/GenomeData.txt");
print $out "ID\tName\tGenes\n";
foreach my $key (keys(%{$data})) {
	print $out $key."\t".$data->{$key}."\t".$genomeHash->{$key}->[0]."\n";
}
close($out);