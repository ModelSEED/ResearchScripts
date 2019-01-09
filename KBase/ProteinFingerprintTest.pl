#!/usr/bin/perl -w

use strict;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $path = $ARGV[0];
my $size = $ARGV[1];
my $location = $ARGV[2];
my $scanrange = $ARGV[3];
my $genomes = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($path."/AllPatricGenomes.txt")}));
my $totalgenomes = keys(%{$genomes});

my $count = 0;
my $hash = {};
my $genomelist = [keys(%{$genomes})];
for (my $i=0; $i < @{$genomelist}; $i++) { 
	my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/vol/patric3/downloads/genomes/".$genomelist->[$i]."/".$genomelist->[$i].".PATRIC.faa");
	my $id = "";
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/>([^\s^\t]+)[\s\t]/) {
			if (defined($id)) {
				my $found = 0;
				my $start = $location - $scanrange;
				if ($start < 0) {
					$start = 0;
				}
				my $stop = $location + $scanrange;
				if ($stop > (length($seq)-1-$size)) {
					$stop = length($seq)-1-$size;
				}
				for (my $k=0; $k < ($stop-$start); $k++) {
					my $query = substr($seq,$start+$k,$size);
					if (defined($hash->{$query})) {
						&add_node($id,$query,$hash,1,$seq,$start+$k,$size);
						$found = 1;
						last;
					}
				}
				if ($found == 0) {
					if (($location+$size) < length($seq)) {
						my $query = substr($seq,$location,$size);
						&add_node($id,$query,$hash,1,$seq,$location,$size);
					}
				}
				$count++;
			}
			$id = $1;
			$seq = "";
		} else {
			$seq .= $lines->[$j];
		}
	}
	if ($count > 1000) {
		last;
	}
}

my $stats = {};
my $output = &compute_stats($stats,$hash,$size,1);
Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/stats.json",[Bio::KBase::ObjectAPI::utilities::TOJSON($stats,1)]);

sub add_node {
	my ($id,$query,$fhash,$level,$seq,$start,$ssize) = @_;
	if ($start+$ssize*($level+1) >= length($seq)) {
		push(@{$fhash->{genes}},$id);
		return 1;
	}
	if ($level <= 14) {
		if (!defined($fhash->{$query})) {
			$fhash->{$query} = {};
		}
		my $newhash = $fhash->{$query};
		my $newquery = substr($seq,$start+$ssize*$level,$ssize);
		return &add_node($id,$newquery,$newhash,$level+1,$seq,$start,$ssize);
	}
	push(@{$fhash->{$query}->{genes}},$id);
	return 1;
}

sub compute_stats {
	my ($fstats,$fhash,$size,$depth) = @_;
	my $list = [keys(%{$fhash})];
	my $output = {
		total => 0,
		genus => {}
	};
	for (my $i=0; $i < @{$list}; $i++) {
		if ($list->[$i] eq "genes") {
			for (my $j=0; $j < @{$fhash->{genes}}; $j++) {
				$output->{total}++;
			}
		} else {
			my $subout = &compute_stats($fstats,$fhash->{$list->[$i]},$size,$depth+1);
			$output->{total} += $subout->{total};
			my $genuslist = [keys(%{$subout->{genus}})];
			for (my $j=0; $j < @{$genuslist}; $j++) {
				if (!defined($output->{genus}->{$genuslist->[$j]})) {
					$output->{genus}->{$genuslist->[$j]} = 0;
				}
				$output->{genus}->{$genuslist->[$j]} += $subout->{genus}->{$genuslist->[$j]};
			}
		}
	}
	my $kmersize = $size*$depth;
	if (!defined($fstats->{$kmersize}->{genedist}->{$output->{total}})) {
		$fstats->{$kmersize}->{genedist}->{$output->{total}} = 0;
	}
	$fstats->{$kmersize}->{genedist}->{$output->{total}}++;
	my $genuscount = keys(%{$output->{genus}});
	if (!defined($fstats->{$kmersize}->{genusdist}->{$genuscount})) {
		$fstats->{$kmersize}->{genusdist}->{$genuscount} = 0;
	}
	$fstats->{$kmersize}->{genusdist}->{$genuscount}++;
	return $output;
}

#Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/hash.json",[Bio::KBase::ObjectAPI::utilities::TOJSON($hash,1)]);