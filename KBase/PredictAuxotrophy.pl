#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Bio::KBase::ObjectAPI::utilities;

my $filearray = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Desktop/ReactionCompoundAssociation.txt");
my $data = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{$filearray}));
my $biocandidate;
foreach my $rxn (keys(%{$data})) {
	my $cpds = [keys(%{$data->{$rxn}})];
	if (@{$cpds} <= 5) {
		foreach my $cpd (@{$cpds}) {
			$biocandidate->{$cpd}->{$rxn} = 1;
		}
	}
}

my $outputfile = "/Users/chenry/Dropbox/workspace/ModelSEED Project/ModelRxnData.json";
my $rxndata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($outputfile)}));
$outputfile = "/Users/chenry/Desktop/BiomassDependencyStatistics.txt";
my $biocpdstats = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($outputfile)}));

my $genomehash;
foreach my $id (keys(%{$rxndata})) {
	foreach my $genome (keys(%{$rxndata->{$id}->{genomes}})) {
		$genomehash->{$genome} = 1;
		my $biodeps = $rxndata->{$id}->{genomes}->{$genome}->{biomassdep};
		foreach my $biodep (@{$biodeps}) {
			$biodep =~ s/_c0$//;
			if (defined($biocandidate->{$biodep}->{$id})) {
				if (!defined($biocpdstats->{$biodep}->{$genome})) {
					$biocpdstats->{$biodep}->{$genome} = {
						gaps => 0,
						all => 0
					};
				}
				$biocpdstats->{$biodep}->{$genome}->{all}++;
				if ($rxndata->{$id}->{genomes}->{$genome}->{gapfill} ne "") {
					$biocpdstats->{$biodep}->{$genome}->{gaps}++;
				}
			}
		}
	}
}

my $resulthash;
foreach my $genome (keys(%{$genomehash})) {
	foreach my $biodep (keys(%{$biocpdstats})) {
		if (defined($biocpdstats->{$biodep}->{$genome})) {
			if ($biocpdstats->{$biodep}->{$genome}->{gaps} > ($biocpdstats->{$biodep}->{gap_ave} + 0.5*$biocpdstats->{$biodep}->{gap_stddev}) || ($biocpdstats->{$biodep}->{$genome}->{gaps} > 0 && $biocpdstats->{$biodep}->{$genome}->{gaps} == $biocpdstats->{$biodep}->{gap_max})) {
				my $stddev = 0;
				if ($biocpdstats->{$biodep}->{gap_stddev} > 0) {
					$stddev = ($biocpdstats->{$biodep}->{$genome}->{gaps} - $biocpdstats->{$biodep}->{gap_ave})/$biocpdstats->{$biodep}->{gap_stddev};
				}
				$resulthash->{$genome}->{$biodep} = {
					stddevs => $stddev,
					gaps => $biocpdstats->{$biodep}->{$genome}->{gaps}
				};
			}
		}	
	}
}

Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Desktop/ReactionCompoundAssociation.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($resulthash)]);
    	