#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use P3DataAPI;
use gjoseqlib;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $path = $ARGV[0];

my $d = P3DataAPI->new();

my $output = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($path."/AllPatricGenomes.txt")}));

$d->{chunk_size} = 5000;

print "Now retrieving rpoB!\n";

my(@fields) = qw(feature_type
	patric_id
	aa_sequence
	product
);
		 
my @res = $d->query("genome_feature",
	["select", @fields],
	["eq", "feature_type", "CDS"],
	["eq", "annotation", "PATRIC"],
	["eq", "product", "DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)"],
);

for (my $i=0; $i < @res; $i++) {
	if ($res[$i]->{patric_id} =~ m/fig\|(\d+\.\d+)\./) {
		my $genomeid = $1;
		if (defined($output->{$genomeid})) {
			push(@{$output->{$genomeid}->{rpoB}},$res[$i]->{aa_sequence});
		}
	}
}

print "Now retrieving rpoC!\n";

@res = $d->query("genome_feature",
	["select", @fields],
	["eq", "feature_type", "CDS"],
	["eq", "annotation", "PATRIC"],
	["eq", "product", "DNA-directed RNA polymerase beta' subunit (EC 2.7.7.6)"],
);

for (my $i=0; $i < @res; $i++) {
	if ($res[$i]->{patric_id} =~ m/fig\|(\d+\.\d+)\./) {
		my $genomeid = $1;
		if (defined($output->{$genomeid})) {
			push(@{$output->{$genomeid}->{rpoB}},$res[$i]->{aa_sequence});
		}
	}
}

Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/AllPatricGenomes2.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($output,1)]);