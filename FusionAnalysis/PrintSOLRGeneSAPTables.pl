#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;
use Data::Dumper;

my $sapsvr = ModelSEED::Client::SAP->new();

my $directory = $ARGV[0];

my $genehash = {};

open(my $fa, "<", $directory."/FusionsTable.txt");
my $fusionhash;
my $fusionheading = <$fa>;
chomp($fusionheading);
while (my $line = <$fa>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	if ($items->[0] ne "Gene") {
		$genehash->{$items->[0]} = [];
		$fusionhash->{$items->[0]} = $line;
	}
}
close($fa);

open(my $fb, "<", $directory."/TrainingTable.txt");
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

my $list = [keys(%{$genehash})];
my $templist = [];
my $headings = ["function","length","location"];
my $genedata;
for (my $i=0; $i < @{$list}; $i++) {
	push(@{$templist},$list->[$i]);
	if (@{$templist} >= 10000) {
		print $i."\n";
		$genedata = $sapsvr->ids_to_data({
			-ids => $templist,
			-data => $headings
		});
		foreach my $id (keys(%{$genedata})) {
			if (defined($genedata->{$id}->[0]->[0])) {
				$genedata->{$id}->[0]->[0] =~ s/\s*#.+//;
				my $array = [split(/\s*;\s+|\s+[\@\/]\s+/,$genedata->{$id}->[0]->[0])];
				$genedata->{$id}->[0]->[0] = join("|",@{$array});
			} else {
				$genedata->{$id}->[0]->[0] = "unknown";
			}
			if (!defined($genedata->{$id}->[0]->[1])) {
				$genedata->{$id}->[0]->[1] = 0;
			}
			if ($genedata->{$id}->[0]->[2] =~ m/(.+)_(\d+)([\+-])(\d+)$/) {
				$genedata->{$id}->[0]->[2] = $1;
				$genedata->{$id}->[0]->[3] = $3;
				my $start = $2;
				my $length = $4;
				my $end = $start+$length;
				if ($genedata->{$id}->[0]->[3] eq "-") {
					$end = $start;
					$start = $end-$length;
				}
				$genedata->{$id}->[0]->[4] = $start;
				$genedata->{$id}->[0]->[5] = $end;
			} else {
				$genedata->{$id}->[0]->[3] = "?"
				$genedata->{$id}->[0]->[4] = "0"
				$genedata->{$id}->[0]->[5] = "0"
			}
			$genehash->{$id} = $genedata->{$id}->[0];
		}
		$genedata = $sapsvr->fids_to_proteins({
			-ids => $templist,
			-sequence => 1
		});
		foreach my $id (keys(%{$genedata})) {
			push(@{$genehash->{$id}},$genedata->{$id});
		}
		$templist = [];
	}
}
$genedata = $sapsvr->ids_to_data({
	-ids => $templist,
	-data => $headings
});
foreach my $id (keys(%{$genedata})) {
	if (defined($genedata->{$id}->[0]->[0])) {
		$genedata->{$id}->[0]->[0] =~ s/\s*#.+//;
		my $array = [split(/\s*;\s+|\s+[\@\/]\s+/,$genedata->{$id}->[0]->[0])];
		$genedata->{$id}->[0]->[0] = join("|",@{$array});
	} else {
		$genedata->{$id}->[0]->[0] = "unknown";
	}
	if (!defined($genedata->{$id}->[0]->[1])) {
		$genedata->{$id}->[0]->[1] = 0;
	}
	if ($genedata->{$id}->[0]->[2] =~ m/(.+)_(\d+)([\+-])(\d+)$/) {
		$genedata->{$id}->[0]->[2] = $1;
		$genedata->{$id}->[0]->[3] = $3;
		my $start = $2;
		my $length = $4;
		my $end = $start+$length;
		if ($genedata->{$id}->[0]->[3] eq "-") {
			$end = $start;
			$start = $end-$length;
		}
		$genedata->{$id}->[0]->[4] = $start;
		$genedata->{$id}->[0]->[5] = $end;
	} else {
		$genedata->{$id}->[0]->[3] = "?"
		$genedata->{$id}->[0]->[4] = "0"
		$genedata->{$id}->[0]->[5] = "0"
	}
	$genehash->{$id} = $genedata->{$id}->[0];
}
$genedata = $sapsvr->fids_to_proteins({
	-ids => $templist,
	-sequence => 1
});
foreach my $id (keys(%{$genedata})) {
	push(@{$genehash->{$id}},$genedata->{$id});
}

open (my $oa, ">", $directory."/SOLR-FusionsTable.txt");
print $oa $fusionheading."\tfunction\tlength\tcontig\tdirection\tstart\tstop\tsequence\n";
foreach my $item (keys(%{$fusionhash})) {
	print $oa $fusionhash->{$item}."\t".join("\t",@{$genehash->{$item}})."\n";
}
close($oa);

open (my $ob, ">", $directory."/SOLR-TrainingTable.txt");
print $ob $trainingheading."\tfunction\tlength\tcontig\tdirection\tstart\tstop\tsequence\n";
foreach my $item (keys(%{$traininghash})) {
	print $ob $traininghash->{$item}."\t".join("\t",@{$genehash->{$item}})."\n";
}
close($ob);