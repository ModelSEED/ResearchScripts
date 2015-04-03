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
	$tmparray->[4] = 0;
	$CDDData->{$tmparray->[0]} = $tmparray;
}
close($fh);

my $CDDSets = {};
my $sethash = {};
open(CDDSETS, "<", $directory."CDDSets.txt");
while (my $line = <CDDSETS>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$CDDSets->{$tmparray->[0]} = {
		count => $tmparray->[1],
		cdds => [split(/;/,$tmparray->[2])]
	};
	for (my $i=0; $i < @{$CDDSets->{$tmparray->[0]}->{cdds}}; $i++) {
		$sethash->{$CDDSets->{$tmparray->[0]}->{cdds}->[$i]} = $tmparray->[0];
	}
}
close(CDDSETS);

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

my $setpairs;
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $genes = {};
	while (my $line = <$fh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		$genes->{$items->[13]}->{$items->[8]} = $items;
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
				my $setids = [];
				if ($cddtwo->[3] > $cddone->[4] || $cddone->[3] > $cddtwo->[4]) {
					$setids->[0] = $sethash->{$cddlist->[$i]};
					$setids->[1] = $sethash->{$cddlist->[$j]};
				}
				if (!defined($setpairs->{$setids->[0]}->{$setids->[1]}->{$gene})) {
					$setpairs->{$setids->[0]}->{$setids->[1]}->{$gene} = 0;
				}
				$setpairs->{$setids->[0]}->{$setids->[1]}->{$gene}++;
				if (!defined($setpairs->{$setids->[1]}->{$setids->[0]}->{$gene})) {
					$setpairs->{$setids->[1]}->{$setids->[0]}->{$gene} = 0;
				}
				$setpairs->{$setids->[1]}->{$setids->[0]}->{$gene}++;
			}
		}
	}
}

open(my $fhh, ">", $directory."CDDSets.txt");
print $fhh "Representative\tNumber core CDDs\tCore CDDs\tSet link count\tLinked sets\n";
foreach my $set (keys(%{$CDDSets})) {
	my $linkcount = 0;
	my $genelink = "";
	if (defined($setpairs->{$set})) {
		$linkcount = keys(%{$setpairs->{$set}});
		foreach my $newset (keys(%{$setpairs->{$set}})) {
			if (length($genelink) > 0) {
				$genelink .= ";";
			}
			my $genecount = keys(%{$setpairs->{$set}->{$newset}}); 
			$genelink .= $newset.":".$genecount;
		}
	}
	print $fhh $set."\t".$CDDSets->{$set}->{count}."\t".join(";",@{$CDDSets->{$set}->{cdds}})."\t".$linkcount."\t".$genelink."\n";
}
close($fhh);

1;