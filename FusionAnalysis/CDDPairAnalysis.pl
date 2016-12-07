#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];

my $cddpairs = {};
open(my $fh, "<", $directory."CDDPairs.txt");
my $line = <$fh>;
while ($line = <$fh>) {
	chomp($line);
	my $tmparray = [split(/\t/,$line)];
	$cddpairs->{$tmparray->[0]}->{$tmparray->[1]} = $tmparray->[2];
}
close($fh);

my $sets = {};
my $inset = {};
my $threshold = 0.5;
foreach my $cddone (keys(%{$cddpairs})) {
	if (!defined($inset->{$cddone})) {
		$inset->{$cddone} = 1;
		$sets->{$cddone} = {$cddone => 1};
		foreach my $cddtwo (keys(%{$cddpairs->{$cddone}})) {
			if (!defined($inset->{$cddtwo})) {
				my $overcount = 0;
				my $totcount = 0;
				foreach my $cddthree (keys(%{$cddpairs->{$cddtwo}})) {
					$totcount += $cddpairs->{$cddtwo}->{$cddthree};
					if (defined($cddpairs->{$cddone}->{$cddthree})) {
						$overcount += $cddpairs->{$cddtwo}->{$cddthree};
					}
				}
				if ($totcount > 0 && $overcount/$totcount > $threshold) {
					$inset->{$cddtwo} = 1;
					$sets->{$cddone}->{$cddtwo} = 1;
				}
			}
		}
		my $continue = 1;
		while ($continue == 1) {
			my $lastcount = keys(%{$sets->{$cddone}});
			my $cddlist = [keys(%{$sets->{$cddone}})];
			foreach my $cddcore (@{$cddlist}) {
				foreach my $cddtwo (keys(%{$cddpairs->{$cddcore}})) {
					if (!defined($inset->{$cddtwo})) {
						my $overcount = 0;
						my $totcount = 0;
						foreach my $cddthree (keys(%{$cddpairs->{$cddtwo}})) {
							$totcount += $cddpairs->{$cddtwo}->{$cddthree};
							if (defined($sets->{$cddone}->{$cddthree})) {
								$overcount += $cddpairs->{$cddtwo}->{$cddthree};
							}
						}
						if ($totcount > 0 && $overcount/$totcount > $threshold) {
							$inset->{$cddtwo} = 1;
							$sets->{$cddone}->{$cddtwo} = 1;
						}
					}
				}
			}
			my $newcount = keys(%{$sets->{$cddone}});
			if ($newcount == $lastcount) {
				$continue = 0;
			}
		}
	}
}

open(my $fhh, ">", $directory."CDDSets.txt");
print $fhh "Representative\tNumber core CDDs\tCore CDDs\n";
foreach my $set (keys(%{$sets})) {
	my $count = keys(%{$sets->{$set}});
	print $fhh $set."\t".$count."\t".join(";",keys(%{$sets->{$set}}))."\n";
}
close($fhh);

1;