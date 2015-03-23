#!/usr/bin/perl

$|=1;
my $genomefile = $ARGV[0];
my $gene = $ARGV[1];

my $fh;
open($fh, "<", $genomefile);
my $line = <$fh>;
my $genes = {};
while (my $line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$genes->{$items->[0]}->{$items->[8]} = $items;
	$genes->{$items->[13]}->{$items->[8]} = $items;
}
close($fh);
if (defined($genes->{$gene})) {
	my $array = [keys(%{$genes->{$gene}})];
	my $sortedarrays = [ sort { $genes->{$gene}->{$a}->[3] <=> $genes->{$gene}->{$b}->[3] } @{$array} ];
	my $sets = [];
	my $inset = {};
	my $last;
	for (my $j=0;$j < @{$sortedarrays}; $j++) {
		if (!defined($inset->{$j})) {
			$last = $j;
			$inset->{$j} = 1;
			my $newset = [$sortedarrays->[$j]];
			for (my $k=$j+1;$k < @{$sortedarrays}; $k++) {
				if ($genes->{$gene}->{$sortedarrays->[$k]}->[3] > $genes->{$gene}->{$sortedarrays->[$last]}->[4]) {
					$last = $k;
					$inset->{$k} = 1;
					push(@{$newset},$sortedarrays->[$k]);
				}
			}
			push(@{$sets},$newset);
		}
	}
	my $genedata = $genes->{$gene}->{$sets->[0]->[0]};
	for (my $j=0; $j < @{$sets}; $j++) {
		my $current = 0;
		for (my $k=0; $k < ($genedata->[1]/3); $k++) {
			if ($current >= @{$sets->[$j]}) {
				print " ";
			} else {
				if ($k < $genes->{$gene}->{$sets->[$j]->[$current]}->[3]) {
					print " ";
				} elsif ($k == $genes->{$gene}->{$sets->[$j]->[$current]}->[3]) {
					print $genes->{$gene}->{$sets->[$j]->[$current]}->[8];
					$k += length($genes->{$gene}->{$sets->[$j]->[$current]}->[8]);
				} elsif ($k <= $genes->{$gene}->{$sets->[$j]->[$current]}->[4]) {
					print "*";
				} else {
					print " ";
					$current++;
				}
			}
		}
		print "\n";
	}
	for (my $k=0; $k < ($genedata->[1]/3); $k++) {
		print "X";	
	}
	print "\n";
}

1;