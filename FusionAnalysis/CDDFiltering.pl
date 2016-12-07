#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];
my $genome = $ARGV[2];

my $array;
open(my $fh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);
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
	close($fh);
	foreach my $gene (keys(%{$genes})) {
		my $genedata = $genes->{$gene};
		my $cddkeys = [keys(%{$genedata})];
		for (my $j=0; $j < @{$cddkeys}; $j++) {
			for (my $k=$j+1; $k < @{$cddkeys}; $k++) {
				my $overstart = $genedata->{$cddkeys->[$j]}->[0];
				if ($genedata->{$cddkeys->[$k]}->[0] > $overstart) {
					$overstart = $genedata->{$cddkeys->[$k]}->[0];
				}
				my $overend = $genedata->{$cddkeys->[$j]}->[1];
				if ($genedata->{$cddkeys->[$k]}->[1] < $overend) {
					$overend = $genedata->{$cddkeys->[$k]}->[1];
				}
				my $overlap = $overend - $overstart;
				if ($overlap > 0) {
					my $overfrac = $overlap/($genedata->{$cddkeys->[$j]}->[1]-$genedata->{$cddkeys->[$j]}->[0]);
					my $overfract = $overlap/($genedata->{$cddkeys->[$k]}->[1]-$genedata->{$cddkeys->[$k]}->[0]);
					if ($overfrac >= 0.75 && $overfract >= 0.75) {
						if ($genedata->{$cddkeys->[$k]}->[5] > $genedata->{$cddkeys->[$j]}->[5]) {
							splice(@{$cddkeys}, $j, 1);
							$j--;
							last;
						} elsif ($genedata->{$cddkeys->[$k]}->[5] < $genedata->{$cddkeys->[$j]}->[5]) {
							splice(@{$cddkeys}, $k, 1);
							$k--;
						}
					}
				}
			}
		}
	}
	open($fhh, ">", $directory.$array->[$i]);
	foreach my $gene (keys(%{$genes})) {
		foreach my $cdd (keys(%{$genes->{$gene}})) {
			print $fhh join("\t",@{$genes->{$gene}->{$cdd}})."\n";
		}
	}
	close($fhh);
}

1;