#!/usr/bin/perl

use strict;

my $directory = $ARGV[0];

my $functionhash;
open(my $fh, "<", $directory."/FunctionAnalysis.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$functionhash->{$items->[0]} = {
		genes => $items->[1],
		genomes => $items->[2],
		fusions => 0,
		fusiongenomes => 0
	};
}
close($fh);

my $Genomes;
my $fusiongenomehash;
open(my $fhh, "<", $directory."/AllFusions.txt");
while (my $line = <$fhh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	my $genomeid;
	if ($items->[0] =~ m/(fig\|\d+\.\d+)\.(.+)/) {
		$genomeid = $1;
		if (!defined($Genomes->{$genomeid})) {
			$Genomes->{$genomeid} = 0;
		}
		$Genomes->{$genomeid}++;
	}
	$items->[2] =~ s/\s*#.+//;
	my $fnarray = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$items->[2]))];
	for (my $j=0; $j < @{$fnarray}; $j++) {
		$fnarray->[$j] = lc($fnarray->[$j]);
		$fnarray->[$j] =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
		$fnarray->[$j] =~ s/\s//g;
		$fnarray->[$j] =~ s/\#.*$//g;
		$fnarray->[$j] =~ s/\(ec\)//g;
		$functionhash->{$fnarray->[$j]}->{fusions}++;
		if (!defined($fusiongenomehash->{$fnarray->[$j]}->{$genomeid})) {
			$fusiongenomehash->{$fnarray->[$j]}->{$genomeid} = 1;
			$functionhash->{$fnarray->[$j]}->{fusiongenomes}++;
		}
	}
}
close($fhh);

open ( my $outh, ">", $directory."/GenomeFusionCount.txt");
foreach my $genome (keys(%{$Genomes})) {
	print $outh $genome."\t".$Genomes->{$genome}."\n";
}
close($outh);

open ( my $out, ">", $directory."/FunctionFusionAnalysis.txt");
print $out "Function\tCount\tGenome count\tFusions\tFusion genomes\n";
foreach my $function (keys(%{$functionhash})) {
	print $out $function."\t".$functionhash->{$function}->{genes}."\t".$functionhash->{$function}->{genomes}."\t".$functionhash->{$function}->{fusions}."\t".$functionhash->{$function}->{fusiongenomes}."\n";
}
close($out);