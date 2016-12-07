#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];
my $window = 5;

my $sapsvr = ModelSEED::Client::SAP->new();

print "Loading blast file list!\n";
my $array;
open(my $fhh, "<", $directory."BlastFileList.txt");
while (my $line = <$fhh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fhh);

print "Loading blast data!\n";
my $genomehash;
my $plantgenehash;
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i." ".$array->[$i]."\n";
	open(my $fhhh, "<", "/homes/seaver/Projects/Plants_PubSEED_Sims/Blast_Results/".$array->[$i]);
	while (my $line = <$fhhh>) {
		chomp($line);
		my $blastdata = [split(/\t/,$line)];
		if ($blastdata->[1] =~ m/(fig\|\d+\.\d+)\./) {
			$genomehash->{$1} = 1;
			$plantgenehash->{$blastdata->[1]}->{$blastdata->[0]} = 1;
		}
	}
	close($fhhh);
}

open(my $f, ">", $directory."GenomeList.txt");
print $f join("\n",keys(%{$genomehash}));
close($f);

open(my $fh, ">", $directory."BlastOutput.txt");
foreach my $gene (keys(%{$plantgenehash})) {
	print $fh $gene."\t".join(";",keys(%{$plantgenehash->{$gene}}))."\n";
}
close($fh);