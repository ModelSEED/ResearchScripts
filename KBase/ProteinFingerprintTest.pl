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
				if ($count > 100) {
					Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/hash.json",[Bio::KBase::ObjectAPI::utilities::TOJSON($hash,1)]);
					exit(0);
				}
				print "ID:".$id."\n";
				print "Sequence:".$seq."\n";
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
						print "0:".$query."\n";
						&add_node($id,$query,$hash,1,$seq,$start+$k,$size);
						$found = 1;
						last;
					}
				}
				if ($found == 0) {
					my $query = substr($seq,$location,$size);
					print "0:".$query."\n";
					&add_node($id,$query,$hash,1,$seq,$location,$size);
				}
				$count++;
			}
			$id = $1;
			$seq = "";
		} else {
			$seq .= $lines->[$j];
		}
	}
}

sub add_node {
	my ($id,$query,$fhash,$level,$seq,$start,$ssize) = @_;
	print $level.":".$query."\n";
	if ($start+$ssize*($level+1) >= length($seq)) {
		print "DONE-EARLY\n";
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
	print "DONE\n";
	push(@{$fhash->{$query}},$id);
	return 1;
}

Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/hash.json",[Bio::KBase::ObjectAPI::utilities::TOJSON($hash,1)]);