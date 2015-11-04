#!/usr/bin/perl

use strict;
use warnings;

my $directory = $ARGV[0];
my $inputdir = $ARGV[1];

my $genehash = {};
my $cddgenehash = {};

open(my $fa, "<", $directory."/SOLR-FusionsTable.txt");
my $fusionhash;
my $fusionheading = <$fa>;
chomp($fusionheading);
while (my $line = <$fa>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$genehash->{$items->[0]} = [];
	$fusionhash->{$items->[0]} = $line;
}
close($fa);

open(my $fb, "<", $directory."/SOLR-TrainingTable.txt");
my $traininghash;
my $trainingheading = <$fb>;
chomp($trainingheading);
while (my $line = <$fb>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$genehash->{$items->[0]} = [];
	$traininghash->{$items->[0]} = $line;
}
close($fb);

open(my $fg, "<", $directory."/SOLR-CDDTable.txt");
my $cddhash = {};
my $cddheading =  <$fg>;
chomp($cddheading);
while (my $line = <$fg>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$cddhash->{$items->[0]} = $line;
}
close($fg);

my $array;
open(my $fc, "<", $inputdir."/GenomeList.txt");
while (my $line = <$fc>) {
	chomp($line);
	push(@{$array},$line);
}
close($fc);
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	open(my $fd, "<", $inputdir."/".$array->[$i]);
	my $line = <$fd>;
	while (my $line = <$fd>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		if (defined($genehash->{$items->[13]})) {
			push(@{$genehash->{$items->[13]}},$items->[2].":".$items->[3].":".$items->[4].":".$items->[5]);
			$cddgenehash->{$items->[2]}->{$items->[13]} = 1;
		}
	}
	close($fd);
}

open (my $oa, ">", $directory."/SOLR-2FusionsTable.txt");
print $oa $fusionheading."\tcdds\n";
foreach my $item (keys(%{$fusionhash})) {
	print $oa $fusionhash->{$item}."\t".join(";",@{$genehash->{$item}})."\n";
}
close($oa);

open (my $ob, ">", $directory."/SOLR-2TrainingTable.txt");
print $ob $trainingheading."\tcdds\n";
foreach my $item (keys(%{$traininghash})) {
	print $oa $traininghash->{$item}."\t".join(";",@{$genehash->{$item}})."\n";
}
close($ob);

open (my $oc, ">", $directory."/SOLR-2CDDTable.txt");
print $oc $cddheading."\tgenes\n";
foreach my $item (keys(%{$cddhash})) {
	print $oc $cddhash->{$item}."\t".join(";",keys(%{$cddgenehash->{$item}}))."\n";
}
close($oc);