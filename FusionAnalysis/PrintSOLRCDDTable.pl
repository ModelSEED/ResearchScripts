#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;

my $directory = $ARGV[0];

my $cddhash = {};
my $cddnames = {};
open(my $fh, "<", $directory."/TableS10-CDDStats.txt");
my $line = <$fh>;
while ($line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$cddnames->{$items->[2]} = $items->[0];
	$cddhash->{$items->[0]} = {
		id => $items->[0],
		accession => undef,
		name => $items->[2],
		"length" => $items->[1],
		genes => $items->[3],
		fullgenes => $items->[4],
		is_full_gene => $items->[5],
		longgenes => $items->[6],
		description => undef,
		set => undef
	};
}
close($fh);

my $svr = Bio::KBase::CDMI::CDMIClient->new_for_script();

my $cdds = $svr->get_entity_ConservedDomainModel(
	[keys(%{$cddhash})],
	["id","accession","short_name","description"]
);
#print Data::Dumper->Dump([$cdds]);
#exit();
foreach my $cdd (keys(%{$cdds})) {
	$cddhash->{$cdd}->{accession} = $cdds->{$cdd}->{accession};
	$cddhash->{$cdd}->{description} = $cdds->{$cdd}->{description};
}

my $cddsethash = {};
open(my $fb, "<", $directory."/CDDSets.txt");
$line = <$fb>;
while ($line = <$fb>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	my $array = [split(/;/,$items->[2])];
	my $finalcdds = [];
	for (my $i=0; $i < @{$array}; $i++) {
		if (length($array->[$i]) > 0) {
			push(@{$finalcdds},$cddnames->{$array->[$i]});
			$cddhash->{$cddnames->{$array->[$i]}}->{set} = $cddnames->{$items->[0]};
		}
	}
	my $linkedarray = [split(/;/,$items->[4])];
	my $finalarray = [];
	for (my $i=0; $i < @{$linkedarray}; $i++) {
		if (length($linkedarray->[$i]) > 0) {
			my $arrayitem = [split(/:/,$linkedarray->[$i])];
			if (length($arrayitem->[0]) > 0) {
				push(@{$finalarray},$cddnames->{$arrayitem->[0]}.":".$arrayitem->[1]);
			}
		}
	}
	$cddhash->{$cddnames->{$items->[0]}}->{set} = $cddnames->{$items->[0]};
	$cddsethash->{$cddnames->{$items->[0]}} = {
		id => $cddnames->{$items->[0]},
		accession => $cddhash->{$cddnames->{$items->[0]}}->{accession},
		name => $cddhash->{$cddnames->{$items->[0]}}->{name},
		"length" => $cddhash->{$cddnames->{$items->[0]}}->{"length"},
		genes => $cddhash->{$cddnames->{$items->[0]}}->{genes},
		fullgenes => $cddhash->{$cddnames->{$items->[0]}}->{fullgenes},
		is_full_gene => $cddhash->{$cddnames->{$items->[0]}}->{is_full_gene},
		longgenes => $cddhash->{$cddnames->{$items->[0]}}->{longgenes},
		description => $cddhash->{$cddnames->{$items->[0]}}->{description},
		cdds => $items->[1],
		cddlist => join(";",@{$finalcdds}),
		links => $items->[3],
		linkedsets => join(";",@{$finalarray})
	};
}
close($fb);

my $headings = [qw(
id
accession
name
length
genes
fullgenes
is_full_gene
longgenes
description
set
)];

open(my $fhh, ">", $directory."/SOLR-CDDTable.txt");
print $fhh join("\t",@{$headings})."\n";
foreach my $cdd (keys(%{$cddhash})) {
	for (my $i=0; $i < @{$headings}; $i++) {
		if ($i > 0) {
			 print $fhh "\t";
		}
		if (defined($cddhash->{$cdd}->{$headings->[$i]})) {
			print $fhh $cddhash->{$cdd}->{$headings->[$i]};
		}
	}
	print $fhh "\n";
}
close($fhh);

$headings = [qw(
id
accession
name
length
genes
fullgenes
is_full_gene
longgenes
description
cdds
cddlist
links
linkedsets
)];

open(my $fhhh, ">", $directory."/SOLR-CDDSetTable.txt");
print $fhhh join("\t",@{$headings})."\n";
foreach my $cdd (keys(%{$cddsethash})) {
	for (my $i=0; $i < @{$headings}; $i++) {
		if ($i > 0) {
			 print $fhhh "\t";
		}
		if (defined($cddsethash->{$cdd}->{$headings->[$i]})) {
			print $fhhh $cddsethash->{$cdd}->{$headings->[$i]};
		}
	}
	print $fhhh "\n";
}
close($fhhh);
