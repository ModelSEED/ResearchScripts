#!/usr/bin/perl

use POSIX;

$|=1;
my $genome = $ARGV[0];
my $gene = $ARGV[1];
my $directory = $ARGV[2];

#Paramters
my $minalignment = 50;
my $mincddfraction = 0.5;
my $minscore = 20;
my $minsglink = 0;

my $CDDData = {};
open(CDD, "<", $directory."CDDList.txt");
while (my $line = <CDD>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$CDDData->{$tmparray->[0]} = $tmparray;
}
close(CDD);

my $CDDSets = {};
my $sethash = {};
open(CDDSETS, "<", $directory."CDDSets.txt");
while (my $line = <CDDSETS>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$CDDSets->{$tmparray->[0]} = {
		count => $tmparray->[1],
		cdds => [split(/;/,$tmparray->[2])],
		linkcount => $tmparray->[3],
		linkedsets => [split(/;/,$tmparray->[4])]
	};
	for (my $i=0; $i < @{$CDDSets->{$tmparray->[0]}->{cdds}}; $i++) {
		$sethash->{$CDDSets->{$tmparray->[0]}->{cdds}->[$i]} = $tmparray->[0];
	}
}
close(CDDSETS);

my $fh;
open($fh, "<", $directory.$genome);
my $line = <$fh>;
my $genes = {};
while (my $line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	if (defined($CDDData->{$items->[8]})) {
		push(@{$genes->{$items->[0]}},$items);
		push(@{$genes->{$items->[13]}},$items);
	}
}

close($fh);
my $maxlength = 200;
if (defined($genes->{$gene})) {
	print "Gene\tLength\tFunction\tDivide\tScore\tLeft\tRight\tOverlap\tLeft SG\tRight SG\tOverlap SG\tMatches\tBest left\tBest right\tBest left align\tBest right align\n";
	#BEGIN SHARED CODE
	my $gd = $genes->{$gene};
	my $startstops = [];
	my $numcdd = @{$gd};
	my $numsg = 0;
	foreach my $cdddata (@{$gd}) {
		$length = $cdddata->[1];
		$function = $cdddata->[6];
		my $cdd = $cdddata->[8];
		my $alignlen = $cdddata->[7];
		my $alignfract = $alignlen/$CDDData->{$cdd}->[2];
		my $score = $cdddata->[5];
		my $sglinks = $CDDData->{$cdd}->[3];
		if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
			$numsg++;
			if (!defined($startstops->[$cdddata->[3]])) {
				$startstops->[$cdddata->[3]] = [0,0];
			}
			if (!defined($startstops->[$cdddata->[4]])) {
				$startstops->[$cdddata->[4]] = [0,0];
			}
			$startstops->[$cdddata->[4]]->[1]++;
			$startstops->[$cdddata->[3]]->[0]++;
		}
	}
	my $bestscore = 0;
	my $position = 0;
	my $leftcdd = 0;
	my $rightcdd = $numsg;
	for (my $i=0; $i < @{$startstops}; $i++) {
		if (defined($startstops->[$i])) {
			$leftcdd += $startstops->[$i]->[1];
			$rightcdd -= $startstops->[$i]->[0];
			if ($leftcdd > 0 && $rightcdd > 0) {
				my $maxleftalign = 0;
				my $maxrightalign = 0;
				foreach my $cdddata (@{$gd}) {
					my $cdd = $cdddata->[8];
					my $alignlen = $cdddata->[7];
					my $alignfract = $alignlen/$CDDData->{$cdd}->[2];
					my $score = $cdddata->[5];
					my $sglinks = $CDDData->{$cdd}->[3];
					if ($cdddata->[4] <= $i) {
						if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
							if ($cdddata->[7] >= $maxleftalign) {
								$maxleftalign = $cdddata->[7];
							}
						}
					} elsif ($cdddata->[3] >= $i) {
						if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
							if ($cdddata->[7] >= $maxrightalign) {
								$maxrightalign = $cdddata->[7];
							}
						}
					}
				}
				if ($maxleftalign+$maxrightalign > $bestscore) {
					$bestscore = $maxleftalign+$maxrightalign;
					$position = $i;
				}
			}
		}
	}
	if ($bestscore > 0) {
		$position++;
		$leftcdd = 0;
		$rightcdd = 0;
		my $leftsg = 0;
		my $rightsg = 0;
		my $numcdd = @{$gd};
		my $matches = 0;
		my $leftcddhash;
		my $rightcddhash;
		my $bestleft = 0;
		my $bestright = 0;
		my $maxleftalign = 0;
		my $maxrightalign = 0;
		foreach my $cdddata (@{$gd}) {
			$length = $cdddata->[1];
			$function = $cdddata->[6];
			my $cdd = $cdddata->[8];
			my $alignlen = $cdddata->[7];
			my $alignfract = $alignlen/$CDDData->{$cdd}->[2];
			my $score = $cdddata->[5];
			my $sglinks = $CDDData->{$cdd}->[3];
			if ($cdddata->[4] <= $position) {
				$leftcdd++;
				if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
					$leftcddhash->{$cdd} = 1;
					$leftsg++;
					if ($CDDData->{$cdd}->[5] >= $bestleft) {
						$bestleft = $CDDData->{$cdd}->[5];
					}
					if ($CDDData->{$cdd}->[7] >= $maxleftalign) {
						$maxleftalign = $CDDData->{$cdd}->[7];
					}
				}
			} elsif ($cdddata->[3] >= $position) {
				$rightcdd++;
				if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
					$rightcddhash->{$cdd} = 1;
					$rightsg++;
					if ($CDDData->{$cdd}->[5] >= $bestright) {
						$bestright = $CDDData->{$cdd}->[5];
					}
					if ($CDDData->{$cdd}->[7] >= $maxrightalign) {
						$maxrightalign = $CDDData->{$cdd}->[7];
					}
				}
			}
		}
		foreach my $cdd (keys(%{$rightcddhash})) {
			if (defined($leftcddhash->{$cdd})) {
				$matches++;
			}
		}
		my $overlap = $numcdd - $leftcdd - $rightcdd;
		my $oversg = $numsg - $leftsg - $rightsg;
		print $gene."\t".$length."\t".$function."\t".$position."\t".$bestscore."\t".$leftcdd."\t".$rightcdd."\t".$overlap."\t".$leftsg."\t".$rightsg."\t".$oversg."\t".$matches."\t".$bestleft."\t".$bestright."\t".$maxleftalign."\t".$maxrightalign."\n";
	}
	#Building sets of nonoverlapping genes to aid in printing of domains
	my $sortedarrays = [ sort { $a->[3] <=> $b->[3] } @{$gd} ];
	my $sets = [];
	my $inset = {};
	my $last;
	for (my $j=0;$j < @{$sortedarrays}; $j++) {
		if (!defined($inset->{$j})) {
			$last = $j;
			$inset->{$j} = 1;
			my $newset = [$sortedarrays->[$j]];
			for (my $k=$j+1;$k < @{$sortedarrays}; $k++) {
				if ($sortedarrays->[$k]->[3] >  $sortedarrays->[$last]->[4]) {
					if (!defined($inset->{$k})) {
						$last = $k;
						$inset->{$k} = 1;
						push(@{$newset},$sortedarrays->[$k]);
					}
				}
			}
			push(@{$sets},$newset);
		}
	}
	my $genedata = $gd->[0];
	my $fraction = $maxlength/($genedata->[1]/3);
	my $first = 1;
	my $setcount = 0;
	my $allcddsets = {};
	for (my $j=0; $j < @{$sets}; $j++) {
		$first = 1;
		my $current = 0;
		for (my $k=0; $k < ($fraction*$genedata->[1]/3); $k++) {
			if ($current >= @{$sets->[$j]}) {
				print " ";
			} else {
				my $currcdd = $sets->[$j]->[$current];
				if ($k < $fraction*$currcdd->[3]) {
					print " ";
				} elsif ($k >= $fraction*$currcdd->[3] && $first == 1) {
					my $score = $currcdd->[5];
					my $alignlen = $currcdd->[7];
					my $cddfract = floor(100*$alignlen/$CDDData->{$currcdd->[8]}->[2]);
					my $sgcount = $CDDData->{$currcdd->[8]}->[3];
					$first = 0;
					my $sg = "N";
					if ($sgcount > $minsglink && $alignlen > $minalignment && $cddfract >= 100*$mincddfraction && $score >= $minscore) {
						$sg = "Y";
					}
					my $cddset = $sethash->{$currcdd->[8]};
					if (!defined($allcddsets->{$cddset})) {
						$allcddsets->{$cddset} = $setcount;
						$setcount++;
					}
					my $linkcount = $CDDSets->{$cddset}->{linkcount};
					$cddset = $allcddsets->{$cddset};
					print $currcdd->[8]."|".$sg.$cddset."(".$score.")".$alignlen."[".$cddfract."]".$sgcount;
					$k += length($currcdd->[8]."|".$sg.$cddset."(".$score.")".$alignlen."[".$cddfract."]".$sgcount);
				} elsif ($k <= $fraction*$currcdd->[4]) {
					print "*";
				} else {
					print " ";
					$first = 1;
					$current++;
				}
			}
		}
		print "\n";
	}
	for (my $k=0; $k < ($fraction*$genedata->[1]/3); $k++) {
		if ($bestscore > 0 && floor($fraction*$position) == $k) {
			print "O";
		} else {
			print "X";
		}
	}
	print "\n";
}

1;