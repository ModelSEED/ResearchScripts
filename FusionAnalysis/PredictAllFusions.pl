#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];

my $array;
open(my $fh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);

$array = ["kb|g.1870"];

for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	system("perl /vol/model-prod/kbase/ResearchScripts/FusionAnalysis/PredictGenomeFusions.pl ".$inputdir." ".$directory." ".$array->[$i]);
}

1;