#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];

my $CDDData = {};
my $cddfunctions = {};
my $array;
open(my $fh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);
my $pairhash = {};
open(my $fhh, "<", $inputdir."CDDPairs.txt");
while (my $line = <$fhh>) {
	chomp($line);
	my $temparray = [split(/\t/,$line)];
	$pairhash->{$temparray->[0]}->{$temparray->[1]} = [];
}
close($fhh);
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $genes = {};
	while (my $line = <$fh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		$genes->{$items->[0]}->{$items->[8]} = $items;
	}
	foreach my $gene (keys(%{$genes})) {
		my $genedata = $genes->{$gene};
		foreach my $cdd (keys(%{$pairhash})) {
			if (defined($genedata->{$cdd})) {
				foreach my $cddd (keys(%{$pairhash->{$cdd}})) {
					my $cdata = $genedata->{$cdd};
					my $cdatat = $genedata->{$cddd};
					if ($cdata->[3] > $cdatat->[4] || $cdatat->[3] > $cdata->[4]) {
						push(@{$pairhash->{$cdd}->{$cddd}},$gene);
					}
				}
			}
		}
	}
}

open(GENELIST, ">", $directory."PairGenes.txt");
print GENELIST "CDD one\tCDD two\tCount\tGenes\n";
foreach my $cdd (keys(%{$pairhash})) {
	foreach my $cddd (keys(%{$pairhash->{$cdd}})) {
		print GENELIST $cdd."\t".$cddd."\t".@{$pairhash->{$cdd}->{$cddd}}."\t".join(";",@{$pairhash->{$cdd}->{$cddd}})."\n";
	}
}
close(GENELIST);

1;