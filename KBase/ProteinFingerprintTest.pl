#!/usr/bin/perl -w

use strict;
#use JSON::XS;
use JSON;
use POSIX;
local $| = 1;

my $path = $ARGV[0];
my $size = $ARGV[1];
my $location = $ARGV[2];
my $scanrange = $ARGV[3];
my $segements = $ARGV[4];
my $genomes = &FROMJSON(join("\n",@{&LOADFILE($path."/AllPatricGenomes.txt")}));
my $totalgenomes = keys(%{$genomes});

my $count = 0;
my $hash = {};
my $funchash = {};
my $genomelist = [keys(%{$genomes})];
for (my $i=0; $i < @{$genomelist}; $i++) { 
	my $lines = &LOADFILE("/vol/patric3/downloads/genomes/".$genomelist->[$i]."/".$genomelist->[$i].".PATRIC.faa");
	my $id;
	my $func;
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/^>([^\s^\t]+)/) {
			my $newid = $1;
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
			if ($lines->[$j] =~ m/^>([^\s^\t]+)\s\s\s(.+)\s\s\s/) {
				$funchash->{$newid} = $2;
			}
			$id = $newid;
			$seq = "";
		} else {
			$seq .= $lines->[$j];
		}
	}
	if ($count > 1000) {
		last;
	}
}
print "Count:".$count."\n";
my $stats = {};
my $output = &compute_stats($stats,$hash,$size,1);
&PRINTFILE($path."/stats.json",[&TOJSON($stats,1)]);

sub add_node {
	my ($id,$query,$fhash,$level,$seq,$start,$ssize) = @_;
	if ($start+$ssize*($level+1) >= length($seq)) {
		push(@{$fhash->{genes}},$id);
		return 1;
	}
	if ($level <= $segements) {
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
		genus => {},
		funcs => {}
	};
	for (my $i=0; $i < @{$list}; $i++) {
		if ($list->[$i] eq "genes") {
			for (my $j=0; $j < @{$fhash->{genes}}; $j++) {
				$output->{total}++;
				if ($fhash->{genes}->[$j] =~ m/\|(\d+\.\d+)\./) {
					my $species = $genomes->{$1}->{n};
					my $array = [split(/\s/,$species)];
					if (!defined($output->{genus}->{$array->[0]})) {
						$output->{genus}->{$array->[0]} = 0;
					}
					$output->{genus}->{$array->[0]}++;
				}
				if (defined($funchash->{$fhash->{genes}->[$j]})) {
					if (!defined($output->{funcs}->{$funchash->{$fhash->{genes}->[$j]}})) {
						$output->{funcs}->{$funchash->{$fhash->{genes}->[$j]}} = 0;
					}
					$output->{funcs}->{$funchash->{$fhash->{genes}->[$j]}}++;
				}
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
			$genuslist = [keys(%{$subout->{funcs}})];
			for (my $j=0; $j < @{$genuslist}; $j++) {
				if (!defined($output->{funcs}->{$genuslist->[$j]})) {
					$output->{funcs}->{$genuslist->[$j]} = 0;
				}
				$output->{funcs}->{$genuslist->[$j]} += $subout->{funcs}->{$genuslist->[$j]};
			}
		}
	}
	my $kmersize = $size*($depth-1);
	if ($kmersize > 0) {
		if (!defined($fstats->{$kmersize}->{genedist}->{$output->{total}})) {
			$fstats->{$kmersize}->{genedist}->{$output->{total}} = 0;
		}
		$fstats->{$kmersize}->{genedist}->{$output->{total}}++;
		my $genuscount = keys(%{$output->{genus}});
		if (!defined($fstats->{$kmersize}->{genusdist}->{$genuscount})) {
			$fstats->{$kmersize}->{genusdist}->{$genuscount} = 0;
		}
		$fstats->{$kmersize}->{genusdist}->{$genuscount}++;
		my $funccount = keys(%{$output->{funcs}});
		if (!defined($fstats->{$kmersize}->{funcdist}->{$funccount})) {
			$fstats->{$kmersize}->{funcdist}->{$funccount} = 0;
		}
		$fstats->{$kmersize}->{funcdist}->{$funccount}++;
		
		my $genuslist = [sort { $output->{genus}->{$a} <=> $output->{genus}->{$b} } keys(%{$output->{genus}})];
		for (my $i=0; $i < @{$genuslist}; $i++) {
			if (defined($genuslist->[$i])) {
				my $fract = $output->{genus}->{$genuslist->[$i]}/$output->{total};
				$fract = floor($fract/10);
				if (!defined($fstats->{$kmersize}->{genusfract}->[$i]->[$fract])) {
					$fstats->{$kmersize}->{genusfract}->[$i]->[$fract] = 0;
				}
				$fstats->{$kmersize}->{genusfract}->[$i]->[$fract]++;
			}
		}
		my $funclist = [sort { $output->{funcs}->{$a} <=> $output->{funcs}->{$b} } keys(%{$output->{funcs}})];
		for (my $i=0; $i < @{$funclist}; $i++) {
			if (defined($funclist->[$i])) {
				my $fract = $output->{funcs}->{$funclist->[$i]}/$output->{total};
				$fract = floor($fract/10);
				if (!defined($fstats->{$kmersize}->{funcfract}->[$i]->[$fract])) {
					$fstats->{$kmersize}->{funcfract}->[$i]->[$fract] = 0;
				}
				$fstats->{$kmersize}->{funcfract}->[$i]->[$fract]++;
			}
		}
	}
	return $output;
}

sub LOADFILE {
    my ($filename) = @_;
    my $DataArrayRef = [];
    open (my $fh, "<", $filename);
    while (my $Line = <$fh>) {
        $Line =~ s/\r//;
        chomp($Line);
        push(@{$DataArrayRef},$Line);
    }
    close($fh);
    return $DataArrayRef;
}

sub PRINTFILE {
    my ($filename,$arrayRef) = @_;
    open ( my $fh, ">", $filename);
    foreach my $Item (@{$arrayRef}) {
    	print $fh $Item."\n";
    }
    close($fh);
}

sub TOJSON {
    my ($ref,$prettyprint) = @_;
    my $JSON = JSON->new->utf8(1);
    if (defined($prettyprint) && $prettyprint == 1) {
		$JSON->pretty(1);
    }
    return $JSON->encode($ref);
}

sub FROMJSON {
    my ($data) = @_;
    return decode_json $data;
}