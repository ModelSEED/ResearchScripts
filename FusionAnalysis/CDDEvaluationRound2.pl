#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];
my $genome = $ARGV[2];
my $CDDData = {};
open(my $fh, "<", $directory."CDDList.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$tmparray->[3] = 0;
	$tmparray->[5] = 0;
	$CDDData->{$tmparray->[0]} = $tmparray;
}
close($fh);
my $cddpairs = {};
my $array;
open(my $gl, "<", $inputdir."GenomeList.txt");
while (my $line = <$gl>) {
	chomp($line);
	push(@{$array},$line);
}
close($gl);
if (defined($genome)) {
	$array = [$genome];
}
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $genes = {};
	while (my $line = <$fh>) {
		chomp($line);
		$genes->{$items->[13]}->{$items->[8]} = $items;
		my $items = [split(/\t/,$line)];
		my $genefraction = 3*$items->[7]/$items->[1];
		my $cddfraction = $items->[7]/$CDDData->{$items->[8]}->[2];
		my $currcdd = $CDDData->{$items->[8]};
		if ($genefraction >= 0.9 && $cddfraction >= 0.9) {
			$currcdd->[3]++;
		}
		if ($genefraction >= 0.9 && ($currcdd->[2]-($items->[1]/3)) >= 50) {
			if ($items->[10] <= 20 || ($currcdd->[2]-$items->[11]) <= 20) {
				$currcdd->[5]++;
			}
		}
	}
	foreach my $gene (keys(%{$genes})) {
		my $genedata = $genes->{$gene};
		my $cddlist = [keys(%{$genedata})];
		for (my $i=0; $i < @{$cddlist}; $i++) {
			my $cddone = $genes->{$gene}->{$cddlist->[$i]};
			for (my $j=$i+1; $j < @{$cddlist}; $j++) {
				my $cddtwo = $genes->{$gene}->{$cddlist->[$j]};
				my $maxstart = $cddone->[3];
				if ($cddtwo->[3] > $maxstart) {
					$maxstart = $cddtwo->[3];
				}
				my $minstop = $cddone->[4];
				if ($cddtwo->[4] < $minstop) {
					$minstop = $cddtwo->[4];
				}
				my $cddonelen = $CDDData->{$cddlist->[$i]}->[2];
				my $cddtwolen = $CDDData->{$cddlist->[$j]}->[2];
				if (($minstop-$maxstart) > 0) {
					if (($minstop-$maxstart)/$cddonelen > 0.5 && ($minstop-$maxstart)/$cddtwolen > 0.5) {
						if (!defined($cddpairs->{$cddlist->[$i]}->{$cddlist->[$j]})) {
							$cddpairs->{$cddlist->[$i]}->{$cddlist->[$j]} = 0;
						}
						$cddpairs->{$cddlist->[$i]}->{$cddlist->[$j]}++;
					}
				}
			}
		}
	}
}

open(CDDLIST, ">", $directory."CDDList.txt");
print CDDLIST "CDD\tMax align\tSize\tSingleGenes\tTotalGenes\tLongGenes\n";
foreach my $key (keys(%{$CDDData})) {
	print CDDLIST $key."\t".$CDDData->{$key}->[1]."\t".$CDDData->{$key}->[2]."\t".$CDDData->{$key}->[3]."\t".$CDDData->{$key}->[4]."\n";
}
close(CDDLIST);

open(CDDPAIRS, ">", $directory."CDDPairs.txt");
print CDDPAIRS "CDDOne\tCDDTwo\tCount\n";
foreach my $cddone (keys(%{$cddpairs})) {
	foreach my $cddtwo (keys(%{$cddpairs->{$cddone}})) {
		print CDDPAIRS $cddone."\t".$cddtwo."\t".$cddpairs->{$cddone}->{$cddtwo}."\n";
	}
}
close(CDDPAIRS);

1;