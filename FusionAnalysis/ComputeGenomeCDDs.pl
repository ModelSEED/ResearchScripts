#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];
my $genome = $ARGV[1];

open(my $fhh, "<", $directory."CDD-Data.txt");
my $line = <$fhh>;
my $cdddata = {};
while ($line = <$fhh>) {
	my $array = [split(/\t/,$line)];
	$cdddata->{$array->[0]} = {
		id => $array->[0],
		len => $array->[1],
		name => $array->[2],
		genes => $array->[3],
		singlegenes => $array->[4],
		longgenes => $array->[5]
	};
}
close($fhh);

(my $obj,my $info) = get_workspace_object($genome);
my $proteins;
my $output = "";
my $ftrs = $obj->{features};
my $genes = {};
for (my $i=0; $i < @{$ftrs}; $i++) {
	if ($ftrs->[$i]->{protein_translation}) {
		$genes->{$ftrs->[$i]->{id}} = $ftrs->[$i];
		$output .= ">".$ftrs->[$i]->{id}."\n".$ftrs->[$i]->{protein_translation}."\n";
	}
}
open(my $fh, ">", $directory.$genome.".fasta");
print $fh $output;
close ($fh);

system("rpsblast -query ".$directory.$genome.".fasta -evalue 0.01 -seg no -outfmt 6 -db /disks/olson/cdd/data/Cdd -out ".$directory.$genome.".cdd");

open(my $fhhhh, ">", $directory.$genome.".cdd.out");
open(my $fhhh, "<", $directory.$genome.".cdd");
print $fhhhh "Gene	Length	CDD	Start	Stop	Identity	Function	Alignlength	CDD name	Protein	CDDStart	CDDStop	E value	SEEDID\n";
while (my $line = <$fhhh>) {
	$array = [split(/\t/,$line)];
	my $itemarray = [split(/\|/,$array->[1])];
	print $fhhhh $array->[0]."\t".$genes->{$array->[0]}->{dna_sequence_length}."\t".$itemarray->[2]."\t".
		$array->[6]."\t".$array->[7]."\t".$array->[2]."\t".$genes->{$array->[0]}->{function}."\t".
		$array->[3]."\t".$cdddata->{$itemarray->[2]}."\t".$genes->{$array->[0]}->{md5}."\t".$array->[8]."\t".$array->[9]."\t".
		$array->[10]."\t".$genes->{$array->[0]}->{id}."\n";
}
close($fhhh);
close($fhhhh);