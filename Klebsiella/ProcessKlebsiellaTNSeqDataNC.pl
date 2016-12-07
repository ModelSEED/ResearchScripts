#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);

my $directory = $ARGV[0];
my $lbthreshold = $ARGV[1];
my $lungthreshold = $ARGV[2];
my $gap = $ARGV[3];
if (!defined($gap)) {
	$gap = 50;
}
my $averageoperonsize = 0;
my $operoncount = 0;
my $lungess = {};
my $lbess = {};
my $lbconflicts = {};
my $lungconflicts = {};
my $genehash;
my $filename = $directory."/TNSeqRawData.tsv";
open(my $fa, "<", $filename);
my $heading = <$fa>;
while (my $line = <$fa>) {
	my $array = [split(/\t/,$line)];
	if ($array->[6] < $lbthreshold) {
		$lbess->{$array->[0]} = 1;
	}
	if ($array->[7] < $lungthreshold) {
		$lungess->{$array->[0]} = 1;
	}
	$genehash->{$array->[0]} = $array;
}
close($fa);

my $order = [sort { $genehash->{$a}->[3] <=> $genehash->{$b}->[3] } keys(%{$genehash})];

my $phenotypes;
my $conflict = 0;
my $lastoperon;
for (my $i=0; $i < @{$order}; $i++) {
	my $key = $order->[$i];
	my $genelist = [$key];
	if ($genehash->{$key}->[5] eq "-") {
		my $count = $i;
		if ($count > 0) {
			my $continue = 1;
			while ($continue == 1) {
				my $end = $genehash->{$order->[$count]}->[3] - $genehash->{$order->[$count]}->[4];
				my $nextstart = $genehash->{$order->[($count-1)]}->[3];
				my $currgap = $end - $nextstart;
				if ($genehash->{$order->[($count-1)]}->[5] eq "+") {
					$continue = 0;
				} elsif ($currgap > $gap) {
					#print "-".$currgap."\n";
					$continue = 0;
				} else {
					push(@{$genelist},$order->[($count-1)]);
				}
				$count--;
			}
		}
	} else {
		my $count = $i;
		if (defined($order->[($count+1)])) {
			my $continue = 1;
			while ($continue == 1) {
				my $end = $genehash->{$order->[$count]}->[3] + $genehash->{$order->[$count]}->[4];
				my $nextstart = $genehash->{$order->[($count+1)]}->[3];
				my $currgap = $nextstart - $end;
				if ($genehash->{$order->[($count+1)]}->[5] eq "-") {
					$continue = 0;
				} elsif ($currgap  > $gap) {
					#print "+".$currgap."\n";
					$continue = 0;
				} else {
					push(@{$genelist},$order->[($count+1)]);
				}
				$count++;
			}
		}
	}
	my $growth = 1;
	my $lunggrowth = 1;
	if ($genehash->{$key}->[6] < $lbthreshold) {
		$growth = 0;
	}
	if ($genehash->{$key}->[7] < $lungthreshold) {
		$lunggrowth = 0;
	}
	for (my $j=1; $j < @{$genelist}; $j++) {
		if (defined($lungess->{$genelist->[$j]})) {
			if ($lunggrowth == 0) {
				$lungconflicts->{$key} = 1;
			}
		}
		if (defined($lbess->{$genelist->[$j]})) {
			if ($growth == 0) {
				$lbconflicts->{$key} = 1;
			}
		}
	}
	push(@{$phenotypes},["ArgonneLBMedia","KBaseMedia",$growth,$key,"none"]);
	push(@{$phenotypes},["LungMedia","chenry:1449636386278",$lunggrowth,$key,"none"]);
	if (@{$genelist} == 1 && defined($lastoperon)) {
		$averageoperonsize += $lastoperon;
		$operoncount++;
	}
	$lastoperon = @{$genelist};
}

$filename = $directory."/UpstreamEssentialityExplanations.tsv";
open(my $fc, ">", $filename);
print $fc "LB upstreams:\t".keys(%{$lbconflicts})."\t".join("\t",keys(%{$lbconflicts}))."\n";
print $fc "Lung upstreams:\t".keys(%{$lungconflicts})."\t".join("\t",keys(%{$lungconflicts}))."\n";
close($fc);

$averageoperonsize += $lastoperon;
$operoncount++;
$averageoperonsize = $averageoperonsize/$operoncount;
print "Operon count:".$operoncount."\n";
print "Average operon size:".$averageoperonsize."\n";
print "LB upstreams:".keys(%{$lbconflicts})."\n";
print "Lung upstreams:".keys(%{$lungconflicts})."\n";
print "LB essentials:".keys(%{$lbess})."\n";
print "Lung essentials:".keys(%{$lungess})."\n";

$filename = $directory."/KlebsiellaPhenotypesSG.tsv";
open(my $fb, ">", $filename);
print $fb "media	mediaws	growth	geneko	addtlCpd\n";
for (my $i=0; $i < @{$phenotypes}; $i++) {
	print $fb join("\t",@{$phenotypes->[$i]})."\n";
}
close($fb);