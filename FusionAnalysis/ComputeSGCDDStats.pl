#!/usr/bin/perl

$|=1;
my $inputdir = $ARGV[0];

my $cddhash = {};
open(my $fhhh, "<", $inputdir."CDDData.txt");
while (my $line = <$fhhh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$cddhash->{$items->[0]} = $items;
}

my $array;
open(my $fhh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fhh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fhh);
my $genecount = 0;
my $genomecount = 0;
my $cddhitcount = 0;
my $nonovercount = 0;
print "Genome\tGene count\tCDD count\tNonover count\n";
for (my $i=0; $i < @{$array}; $i++) {
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $currcddcount = 0;
	my $currnonovercount = 0;
	my $genehash = {};
	while (my $line = <$fh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		if (defined($cddhash->{$items->[2]}) && $cddhash->{$items->[2]}->[4] > 0 && $items->[7]/$cddhash->{$items->[2]}->[1] >= 0.5) {
			if (!defined($genehash->{$items->[0]})) {
				$genehash->{$items->[0]} = {
					firststop => $items->[4],
					laststart => $items->[3]
				};
			} else {
				if ($genehash->{$items->[0]}->{firststop} > $items->[4]) {
					$genehash->{$items->[0]}->{firststop} = $items->[4];
				}
				if ($genehash->{$items->[0]}->{laststart} < $items->[3]) {
					$genehash->{$items->[0]}->{laststart} = $items->[3];
				}
				if ($genehash->{$items->[0]}->{firststop} < $genehash->{$items->[0]}->{laststart} && !defined($genehash->{$items->[0]}->{nonover})) {
					$genehash->{$items->[0]}->{nonover} = 1;
					$currnonovercount++;
				}
			}
			$currcddcount++;
		}
	}	
	if ($currcddcount > 0) {
		$genomecount++;
	}
	my $currentgenecount = keys(%{$genehash});
	print $array->[$i]."\t".$currentgenecount."\t".$currcddcount."\t".$currnonovercount."\n";
	$genecount += $currentgenecount;
	$cddhitcount += $currcddcount;
	$nonovercount += $currnonovercount;
	close($fh);
}
print "Genomes\t".$genomecount."\nGenes\t".$genecount."\nCDD hits:".$cddhitcount."\nNonover count:".$nonovercount."\n";

1;