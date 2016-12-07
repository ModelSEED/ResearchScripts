#!/usr/bin/perl

$|=1;
my $input = $ARGV[0];
my $output = $ARGV[1];
my $genome = $ARGV[2];

#Paramters
my $minalignment = 50;#50
my $mincddfraction = 0.5;#0.6
my $minscore = 20;#25
my $minsglink = 0;#0

my $CDDData;
my $SingleGeneCDDs = {};
my $LongCDDs = {};
my $genedata = {};
open(my $fh, "<", $output."CDDList.txt");
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
open(CDDSETS, "<", $output."CDDSets.txt");
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

open(my $fhh, "<", $input.$genome);
my $line = <$fhh>;
my $GeneCDDs = {};
my $newcdds = {};
while (my $line = <$fhh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	if (defined($CDDData->{$items->[8]})) {
		push(@{$GeneCDDs->{$items->[13]}},$items);
	}
}

open($fh, ">> ".$output."AllFusions.txt"); 
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
			my $reason = "";
			my $predictedfusion = 0;
			my $leftdist = $length/3-$position;
			if ($length >= 600) {*
				if ($leftcdd > 0) {*
					if ($rightcdd > 0) {*
						if ($matches <= 1000) {*
							if ($maxleftalign > 50) {*
								if ($maxrightalign > 50) {*
									if ($bestleft >= 18) {
										if ($bestright >= 18) {
											if ((3*($maxleftalign+$maxrightalign)/$length) >=0.4) {*
												if ($oversg >= 0 && $position >= 58) {*
													if ($leftdist >= 60) {*
														if (($oversg/($rightsg+$leftsg)) <= 1) {*
															if (($leftlinks+$rightlinks) < 3000) {*
																if ((3*$position/$length) >= 0.1) {*
																	if ((3*$leftdist/$length) >= 0.1) {*
																		if ($count > 1) {*
																			if ($rightlinks < 1500 ) {*
																				if ($leftlinks < 1500) {*
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
				print $fh $gene."\t".$length."\t".$function."\t".$position."\t".$bestscore."\t".$leftcdd."\t".$rightcdd."\t".$overlap."\t".$leftsg."\t".$rightsg."\t".$oversg."\t".$matches."\t".$bestleft."\t".$bestright."\t".$maxleftalign."\t".$maxrightalign."\t".$leftlinks."\t".$rightlinks."\t".$count."\n";
			}
		}
	}
}
close($fh);

1;