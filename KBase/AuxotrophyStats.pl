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

my $biocpdhash = {};
foreach my $id (keys(%{$rxndata})) {
	foreach my $genome (keys(%{$rxndata->{$id}->{genomes}})) {
		my $biodeps = $rxndata->{$id}->{genomes}->{$genome}->{biomassdep};
		foreach my $biodep (@{$biodeps}) {
			$biodep =~ s/_c0$//;
			if (defined($biocandidate->{$biodep}->{$id})) {
				if (!defined($biocpdhash->{$biodep}->{$genome})) {
					$biocpdhash->{$biodep}->{$genome} = {
						gaps => 0,
						all => 0
					};
				}
				$biocpdhash->{$biodep}->{$genome}->{all}++;
				if ($rxndata->{$id}->{genomes}->{$genome}->{gapfill} ne "") {
					$biocpdhash->{$biodep}->{$genome}->{gaps}++;
				}
			}
		}
	}
}

my $resulthash = {};
foreach my $biocpd (keys(%{$biocpdhash})) {
	$resulthash->{$biocpd} = {
		all_ave => 0,
		gap_ave => 0,
		all_max => undef,
		gap_max => undef,
		all_min => undef,
		gap_min => undef,
		all_stddev => 0,
		gap_stddev => 0,
		all_count => 0,
		gap_count => 0,
	};
	foreach my $genome (keys(%{$biocpdhash->{$biocpd}})) {
		$resulthash->{$biocpd}->{all_count}++;
		$resulthash->{$biocpd}->{gap_count}++;
		$resulthash->{$biocpd}->{all_ave} += $biocpdhash->{$biocpd}->{$genome}->{all};
		$resulthash->{$biocpd}->{gap_ave} += $biocpdhash->{$biocpd}->{$genome}->{gaps};
		if (!defined($resulthash->{$biocpd}->{all_max}) || $resulthash->{$biocpd}->{all_max} < $biocpdhash->{$biocpd}->{$genome}->{all}) {
			$resulthash->{$biocpd}->{all_max} = $biocpdhash->{$biocpd}->{$genome}->{all};
		}
		if (!defined($resulthash->{$biocpd}->{gap_max}) || $resulthash->{$biocpd}->{gap_max} < $biocpdhash->{$biocpd}->{$genome}->{gaps}) {
			$resulthash->{$biocpd}->{gap_max} = $biocpdhash->{$biocpd}->{$genome}->{gaps};
		}
		if (!defined($resulthash->{$biocpd}->{all_min}) || $resulthash->{$biocpd}->{all_min} > $biocpdhash->{$biocpd}->{$genome}->{all}) {
			$resulthash->{$biocpd}->{all_min} = $biocpdhash->{$biocpd}->{$genome}->{all};
		}
		if (!defined($resulthash->{$biocpd}->{gap_min}) || $resulthash->{$biocpd}->{gap_min} > $biocpdhash->{$biocpd}->{$genome}->{gaps}) {
			$resulthash->{$biocpd}->{gap_min} = $biocpdhash->{$biocpd}->{$genome}->{gaps};
		}
	}
	$resulthash->{$biocpd}->{all_ave} = $resulthash->{$biocpd}->{all_ave}/$resulthash->{$biocpd}->{all_count};
	$resulthash->{$biocpd}->{gap_ave} = $resulthash->{$biocpd}->{gap_ave}/$resulthash->{$biocpd}->{gap_count};
	foreach my $genome (keys(%{$biocpdhash->{$biocpd}})) {
		$resulthash->{$biocpd}->{all_stddev} += ($biocpdhash->{$biocpd}->{$genome}->{all}-$resulthash->{$biocpd}->{all_ave})*($biocpdhash->{$biocpd}->{$genome}->{all}-$resulthash->{$biocpd}->{all_ave});
		$resulthash->{$biocpd}->{gap_stddev} += ($biocpdhash->{$biocpd}->{$genome}->{gaps}-$resulthash->{$biocpd}->{gap_ave})*($biocpdhash->{$biocpd}->{$genome}->{gaps}-$resulthash->{$biocpd}->{gap_ave});
	}
	$resulthash->{$biocpd}->{all_stddev} = $resulthash->{$biocpd}->{all_stddev}/$resulthash->{$biocpd}->{all_count};
	$resulthash->{$biocpd}->{all_stddev} = sqrt($resulthash->{$biocpd}->{all_stddev});
	$resulthash->{$biocpd}->{gap_stddev} = $resulthash->{$biocpd}->{gap_stddev}/$resulthash->{$biocpd}->{gap_count};
	$resulthash->{$biocpd}->{gap_stddev} = sqrt($resulthash->{$biocpd}->{gap_stddev});
}

Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Desktop/BiomassDependencyStatistics.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($resulthash)]);
    	