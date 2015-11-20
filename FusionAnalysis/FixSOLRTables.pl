#!/usr/bin/perl

use strict;
use warnings;

my $directory = $ARGV[0];

open(my $fc, "<", $directory."/SOLR-CDDTable.txt");
my $line = <$fc>;
chomp($line);
my $heading = $line;
my $cddlist = [];
my $cddhash = {};
my $idlist = [];
while ($line = <$fc>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	push(@{$idlist},$items->[0]);
	$cddhash->{$items->[0]} = $items->[2];
	push(@{$cddlist},$line);
}
close($fc);

open(my $fd, ">", $directory."/SOLR-CDDTable2.txt");
print $fd $heading."\tset_name\n";
for (my $i=0; $i < @{$cddlist}; $i++) {
	print $fd $cddlist->[$i]."\t".$cddhash->{$idlist->[$i]}."\n";
}
close($fd);

open(my $fa, "<", $directory."/SOLR-FusionsTable.txt");
$line = <$fa>;
chomp($line);
my @items = [split(/\t/,$line)];
@items = splice(@items,2,1);
@items = splice(@items,1,1);
open(my $fb, ">", $directory."/SOLR-FusionsTable2.txt");
print $fb join("\t",@items)."\n";
while ($line = <$fa>) {
	chomp($line);
	@items = split(/\t/,$line);
	my $columns = @items;
	if ($columns == 28) {
		@items = splice(@items,2,1);
		@items = splice(@items,1,1);
		my $cdds = [split(/;/,$items[25])];
		for (my $i=0; $i < @{$cdds}; $i++) {
			my $cdddata = [split(/:/,$cdds->[$i])];
			unshift(@{$cdddata},$cddhash->{$cdddata->[0]});
			$cdds->[$i] = join(":",@{$cdddata});
		}
		$items[25] = join(";",@{$cdds});
		for (my $i=0;$i < @items; $i++) {
			if (!defined($items[$i])) {
				$items[$i] = "";
			}
		}
		print $fb join("\t",@items)."\n";
	} else {
		print $line."\n";
	}
}
close($fa);
close($fb);

open(my $fe, "<", $directory."/SOLR-CDDSets.txt");
$line = <$fe>;
open(my $ff, ">", $directory."/SOLR-CDDSets2.txt");
print $ff $line;
while ($line = <$fe>) {
	chomp($line);
	@items = split(/\t/,$line);
	my $columns = @items;
	if ($columns == 13) {
		my $cdds = [split(/;/,$items[10])];
		for (my $i=0; $i < @{$cdds}; $i++) {
			my $cdddata = [split(/:/,$cdds->[$i])];
			unshift(@{$cdddata},$cddhash->{$cdddata->[0]});
			$cdds->[$i] = join(":",@{$cdddata});
		}
		$items[10] = join(";",@{$cdds});
		print $ff join("\t",@items)."\n";
	} else {
		print $line."\n";
	}
}
close($fe);
close($ff);