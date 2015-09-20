#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];

#Paramters
my $minalignment = 50;#50
my $mincddfraction = 0.5;#0.6
my $minscore = 20;#25
my $minsglink = 0;#0

my $FusionList = [qw(
b0002
b0067
b0149
b0153
b0414
b0448
b0449
b0465
b0679
b0695
b0731
b0760
b0794
b0820
b0829
b0879
b0886
b0887
b0914
b0924
b0949
b0993
b1014
b1084
b1101
b1198
b1241
b1262
b1263
b1378
b1387
b1439
b1496
b1513
b1621
b1629
b1694
b1755
b1816
b1817
b1883
b1888
b1900
b2022
b2026
b2063
b2149
b2167
b2169
b2180
b2211
b2213
b2218
b2220
b2255
b2324
b2341
b2370
b2383
b2463
b2468
b2519
b2547
b2552
b2554
b2584
b2599
b2600
b2786
b2818
b2829
b2836
b2881
b2887
b2988
b3052
b3053
b3056
b3072
b3168
b3210
b3352
b3368
b3396
b3409
b3418
b3464
b3486
b3567
b3599
b3639
b3722
b3730
b3746
b3749
b3846
b3863
b3868
b3899
b3940
b3947
b4004
b4006
b4035
b4087
b4136
b4340
b4391
b0025
b0054
b0066
b0156
b0192
b0207
b0529
b0550
b0578
b0600
b0657
b0759
b0808
b0812
b0948
b0980
b1814
b1907
b2049
b2194
b2249
b2286
b2507
b2508
b2557
b2797
b2845
b2875
b2972
b3012
b4471
b3114
b3228
b3335
b3551
b3615
b3765
b3942
b3973
b4090
b4159
b4167
b4390
fig|382464.3.peg.4351
fig|326426.4.peg.1408
fig|644966.3.peg.1460
fig|608534.3.peg.920
fig|555079.3.peg.2426
fig|481448.7.peg.1889
fig|390333.6.peg.218
fig|272557.1.peg.1396
fig|178306.1.peg.2190
fig|266779.9.peg.167
fig|138119.3.peg.3934
fig|478801.5.peg.881
fig|373903.5.peg.523
fig|180281.4.peg.761
fig|226186.1.peg.1442
fig|3702.11.peg.8629
fig|234826.3.peg.438
fig|219305.4.peg.3726
fig|411465.10.peg.1479
fig|342108.5.peg.2398
fig|269800.4.peg.2728
fig|452863.6.peg.2446
fig|103690.1.peg.3826
fig|64091.1.peg.1812
fig|357808.3.peg.923
fig|3702.11.peg.4649
fig|546274.4.peg.1959
fig|4932.3.peg.877
fig|4896.1.peg.3261
fig|717606.6.peg.1174
fig|298653.4.peg.1987
fig|272559.3.peg.965
fig|590998.5.peg.1658
fig|395965.4.peg.2584
fig|469616.3.peg.1833
fig|3702.11.peg.4621
fig|330779.3.peg.297
fig|233413.1.peg.32
fig|267747.1.peg.1824
fig|290400.10.peg.2898
fig|321955.3.peg.1718
fig|391625.5.peg.3242
fig|7227.3.peg.13624
fig|178306.1.peg.2192
fig|243232.1.peg.949
fig|3702.11.peg.9804
fig|381764.6.peg.1714
fig|177439.1.peg.1605
fig|187303.17.peg.2469
fig|521095.6.peg.109
fig|63737.4.peg.3651
fig|313606.3.peg.7087
fig|3702.11.peg.16955
fig|96561.3.peg.1759
fig|391603.3.peg.524
fig|395495.3.peg.2329
fig|272634.1.peg.300
fig|500635.8.peg.2116
fig|4896.1.peg.4275
fig|3702.11.peg.466
fig|349741.3.peg.1549
fig|233150.3.peg.42
fig|376686.6.peg.1713
fig|554065.3.peg.6386
fig|1148.1.peg.2723
fig|316067.3.peg.130
fig|591167.6.peg.6429
fig|3702.11.peg.22518
fig|866775.3.peg.61
fig|572480.3.peg.1028
fig|523794.5.peg.1959
fig|10090.3.peg.25319
fig|69014.3.peg.309
fig|471875.6.peg.815
fig|272635.1.peg.335
fig|362976.10.peg.504
fig|411474.6.peg.1072
fig|498211.3.peg.425
fig|431947.6.peg.153
fig|161528.3.peg.1164
fig|59931.3.peg.1325
fig|83333.1.peg.25
fig|240016.6.peg.2924
fig|177437.4.peg.854
fig|563038.3.peg.1439
fig|521095.6.peg.839
fig|64091.1.peg.2192
fig|243277.1.peg.62
fig|592015.5.peg.1366
fig|479436.6.peg.307
fig|632772.3.peg.4758
fig|479435.6.peg.314
fig|585501.3.peg.1201
fig|528347.5.peg.721
fig|83333.1.peg.4298
fig|330779.3.peg.1576
fig|1127122.3.peg.1544
fig|273057.1.peg.2
fig|438753.3.peg.3255
fig|100226.1.peg.4062
fig|338963.3.peg.1505
fig|525246.3.peg.1726
fig|264198.3.peg.721
fig|314265.3.peg.941
fig|257309.1.peg.1719
fig|639282.3.peg.2135
fig|3702.11.peg.14378
fig|251221.1.peg.2647
fig|272568.11.peg.1069
fig|3702.7.peg.25891
fig|378806.7.peg.2871
fig|469382.4.peg.3193
fig|471853.5.peg.717
fig|515620.4.peg.2275
fig|309798.3.peg.1142
fig|439235.3.peg.5130
fig|367737.4.peg.1188
fig|765698.3.peg.4766
fig|708616.3.peg.127
fig|452638.3.peg.1662
fig|272561.1.peg.628
fig|431947.6.peg.1037
fig|264201.15.peg.1181
fig|75379.4.peg.581
fig|247156.1.peg.747
fig|4932.3.peg.5208
fig|668336.4.peg.1151
fig|3702.11.peg.2993
fig|525146.3.peg.1509
fig|1001739.3.peg.2350
fig|240015.3.peg.2233
fig|471854.4.peg.5604
fig|4932.3.peg.5667
fig|471853.5.peg.728
fig|500633.7.peg.413
fig|188626.3.peg.1493
fig|223926.6.peg.2972
fig|292805.3.peg.329
fig|101510.15.peg.5449
fig|196164.1.peg.1591
fig|1148.1.peg.1732
fig|203267.1.peg.586
fig|519442.4.peg.1288
fig|546271.3.peg.778
fig|405948.11.peg.5492
fig|469381.4.peg.1963
fig|426368.9.peg.603
fig|391774.5.peg.1943
fig|9606.3.peg.25573
)];
my $fusionhash = {};
for (my $i=0; $i < @{$FusionList}; $i++) {
	$fusionhash->{$FusionList->[$i]} = 0;
}

my $filteredfusions = 0;
my $unfilteredfusions = 0;
my $filterfusions = 0;
my $missedfusions = 0;
my $fpfusions = 0;
my $correctfusions = 0;

my $CDDData;
my $SingleGeneCDDs = {};
my $LongCDDs = {};
my $genedata = {};
open(my $fh, "<", $directory."CDDList.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$CDDData->{$tmparray->[0]} = $tmparray;
	if ($tmparray->[5] >= 0) {
		$LongCDDs->{$tmparray->[0]} = 1;
	}
	if ($tmparray->[3] >= 0) {
		$SingleGeneCDDs->{$tmparray->[0]} = 1;
	}
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
		cdds => [split(/;/,$tmparray->[2])],
		linkcount => $tmparray->[3],
		linkedsets => 
	};
	my $array = [split(/;/,$tmparray->[4])];
	my $total = 1;
	for (my $i=0; $i < @{$array}; $i++) {
		my $arraytwo = [split(/:/,$array->[$i])];
		$CDDSets->{$tmparray->[0]}->{$arraytwo->[0]} = [$arraytwo->[1],$arraytwo->[1]];
		$total += $arraytwo->[1];
	}
	foreach my $cddset (keys(%{$CDDSets->{$tmparray->[0]}})) {
		$CDDSets->{$tmparray->[0]}->{$cddset}->[1] = $CDDSets->{$tmparray->[0]}->{$cddset}->[1]/$total;
	}
	for (my $i=0; $i < @{$CDDSets->{$tmparray->[0]}->{cdds}}; $i++) {
		$sethash->{$CDDSets->{$tmparray->[0]}->{cdds}->[$i]} = $tmparray->[0];
	}
}
close(CDDSETS);

open(my $fhh, "<", $directory."GeneCDDHits.txt");
my $line = <$fhh>;
my $GeneCDDs = {};
my $newcdds = {};
while (my $line = <$fhh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	if (defined($CDDData->{$items->[8]})) {
		push(@{$GeneCDDs->{$items->[13]}},$items);
	} elsif (defined($newcdds->{$items->[8]})) {
		my $fraction = 3*$items->[7]/$items->[1];
		if ($items->[7] > $newcdds->{$items->[8]}->[2]) {
			$newcdds->{$items->[8]}->[2] = $items->[7];
		}
		if ($items->[7] > $newcdds->{$items->[8]}->[1]) {
			$newcdds->{$items->[8]}->[1] = $items->[7];
		}	
		if ($fraction > 0.85) {
			$newcdds->{$items->[8]}->[3]++;
		}
		$newcdds->{$items->[8]}->[4]++;
	} else {
		my $fraction = 3*$items->[7]/$items->[1];
		$newcdds->{$items->[8]} = [$items->[8],$items->[7],$items->[7],0,1,0];
		if ($fraction > 0.85) {
			$newcdds->{$items->[8]}->[3]++;
		}
	}	
}
foreach my $cdd (keys(%{$newcdds})) {
	#print join("\t",@{$newcdds->{$cdd}})."\n";
}

open($fh, "> ".$directory."AllFusions.txt"); 
print $fh "Gene\tLength\tFunction\tDivide\tScore\tLeft\tRight\tOverlap\tLeft SG\tRight SG\tOverlap SG\tMatches\tBest left\tBest right\tBest left align\tBest right align\tLeft links\tRight links\tSet count\n";
foreach my $gene (keys(%{$GeneCDDs})) {
	my $gd = $GeneCDDs->{$gene};
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
		my $leftlinks = 0;
		my $rightlinks = 0;
		my $genesethash;
		foreach my $cdddata (@{$gd}) {
			$length = $cdddata->[1];
			$function = $cdddata->[6];
			my $cdd = $cdddata->[8];
			$genesethash->{$sethash->{$cdd}} = 1;
			my $alignlen = $cdddata->[7];
			my $alignfract = $alignlen/$CDDData->{$cdd}->[2];
			my $score = $cdddata->[5];
			my $sglinks = $CDDData->{$cdd}->[3];
			if ($cdddata->[4] <= $position) {
				$leftcdd++;
				if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
					if ($leftlinks < $CDDSets->{$sethash->{$cdd}}->{linkcount}) {
						$leftlinks = $CDDSets->{$sethash->{$cdd}}->{linkcount};
					}
					$leftcddhash->{$cdd} = 1;
					$leftsg++;
					if ($cdddata->[5] >= $bestleft) {
						$bestleft = $cdddata->[5];
					}
					if ($alignlen >= $maxleftalign) {
						$maxleftalign = $alignlen;
					}
				}
			} elsif ($cdddata->[3] >= $position) {
				$rightcdd++;
				if ($sglinks > $minsglink && $alignlen > $minalignment && $alignfract >= $mincddfraction && $score >= $minscore) {
					if ($rightlinks < $CDDSets->{$sethash->{$cdd}}->{linkcount}) {
						$rightlinks = $CDDSets->{$sethash->{$cdd}}->{linkcount};
					}
					$rightcddhash->{$cdd} = 1;
					$rightsg++;
					if ($score >= $bestright) {
						$bestright = $score;
					}
					if ($alignlen >= $maxrightalign) {
						$maxrightalign = $alignlen;
					}
				}
			}
		}
		foreach my $cdd (keys(%{$rightcddhash})) {
			if (defined($leftcddhash->{$cdd})) {
				$matches++;
			}
		}
		my $count = keys(%{$genesethash});
		my $overlap = $numcdd - $leftcdd - $rightcdd;
		my $oversg = $numsg - $leftsg - $rightsg;
		if ($leftcdd > 0 && $rightcdd > 0) {
			if (defined($fusionhash->{$gene})) {
				$fusionhash->{$gene} = 1;
				$unfilteredfusions++;
			}
			my $reason = "";
			my $predictedfusion = 0;
			my $leftdist = $length/3-$position;
			if ($length >= 600) {
				if ($leftcdd > 0) {
					if ($rightcdd > 0) {
						if ($matches <= 1000) {
							if ($maxleftalign > 50) {
								if ($maxrightalign > 50) {
									if ($bestleft >= 18) {
										if ($bestright >= 18) {
											if ((3*($maxleftalign+$maxrightalign)/$length) >=0.4) {
												if ($oversg >= 0 && $position >= 58) {
													if ($leftdist >= 60) {
														if (($oversg/($rightsg+$leftsg)) <= 1) {
															if (($leftlinks+$rightlinks) < 3000) {
																if ((3*$position/$length) >= 0.1) {
																	if ((3*$leftdist/$length) >= 0.1) {
																		if ($count > 1) {
																			if ($rightlinks < 1500 ) {
																				if ($leftlinks < 1500) {
																					$predictedfusion = 1;
																				} else {
																					$reason = "Left links > 1500";
																				}
																			} else {
																				$reason = "Right links > 1500";
																			}
																		} else {
																			$reason = "Count < 1";
																		}
																	} else {
																		$reason = "Right dist too short";
																	}
																} else {
																	$reason = "Left dist too short";
																}
															} else {
																$reason = "Sum of links > 3000";
															}
														} else {
															$reason = "Too many overlapping SG";
														}
													} else {
														$reason = "Position < 60 from right";
													}
												} else {
													$reason = "Position < 60 from left";
												}
											} else {
												$reason = "Total align < 0.4";
											}
										} else {
											$reason = "Best right less than 25";
										}
									} else {
										$reason = "Best left less than 25";
									}
								} else {
									$reason = "Max right less than 50";
								}
							} else {
								$reason = "Max left less than 50";
							}
						} else {
							$reason = "More than 162 matches";
						}
					} else {
						$reason = "Max right less than 50";
					}
				} else {
					$reason = "No right CDD";
				}
			} else {
				$reason = "No left CDD";
			}
			if ($predictedfusion == 1) {
				if (defined($fusionhash->{$gene})) {
					$correctfusions++;
					$fusionhash->{$gene} = 2;
				} else {
					$fpfusions++;
				}
			}
			print $fh $gene."\t".$length."\t".$function."\t".$position."\t".$bestscore."\t".$leftcdd."\t".$rightcdd."\t".$overlap."\t".$leftsg."\t".$rightsg."\t".$oversg."\t".$matches."\t".$bestleft."\t".$bestright."\t".$maxleftalign."\t".$maxrightalign."\t".$leftlinks."\t".$rightlinks."\t".$count."\t".$predictedfusion."\t".$reason."\n";
		}
	}
}
close($fh);
my $missedlist = [];
for (my $i=0; $i < @{$FusionList}; $i++) {
	if ($fusionhash->{$FusionList->[$i]} == 0) {
		$filteredfusions++;
		push(@{$missedlist},$FusionList->[$i]);
	} elsif ($fusionhash->{$FusionList->[$i]} == 1) {
		$missedfusions++;
	}
}
print "Missed:".join("\n",@{$missedlist})."\n";
print "Filtered fusions:".$filteredfusions."\n";
print "Unfiltered fusions:".$unfilteredfusions."\n";
print "Correct fusions:".$correctfusions."\n";
print "Missed fusions:".$missedfusions."\n";
print "False positives:".$fpfusions."\n";

1;