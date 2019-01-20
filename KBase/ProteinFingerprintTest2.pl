#!/usr/bin/perl -w

use strict;
#use JSON::XS;
use JSON;
use POSIX;
use Digest::MD5;
local $| = 1;

my $path = $ARGV[0];
my $size = $ARGV[1];
my $genomes = &FROMJSON(join("\n",@{&LOADFILE($path."/AllPatricGenomes.txt")}));
my $totalgenomes = keys(%{$genomes});

my $count = 0;
my $tooshort = 0;
my $overlapping_sketches = {
	"000" => 0,
	"100" => 0,
	"010" => 0,
	"001" => 0,
	"110" => 0,
	"011" => 0,
	"101" => 0,
	"111" => 0
};
my $hash = {};
my $funchash = {};
my $genomelist = [keys(%{$genomes})];
print "Sketch start:".time()."\n";
for (my $i=0; $i < @{$genomelist}; $i++) { 
	print "Sketching:".$i." of ".$totalgenomes."\n";
	my $lines = &LOADFILE("/vol/patric3/downloads/genomes/".$genomelist->[$i]."/".$genomelist->[$i].".PATRIC.faa");
	my $id;
	my $func;
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/^>([^\s^\t]+)/) {
			my $newid = $1;
			if (defined($id)) {
				if (50+$size < length($seq)) {
					my $overlap = "";
					my $sketch = substr($seq,49,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap = "1";
					} else {
						$overlap = "0";
					}
					push(@{$hash->{$sketch}},[$id,0]);
					$sketch = substr($seq,length($seq)-49-$size,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap .= "1";
					} else {
						$overlap .= "0";
					}
					push(@{$hash->{$sketch}},[$id,1]);
					$sketch = substr($seq,length($seq)/2-$size/2,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap .= "1";
					} else {
						$overlap .= "0";
					}
					$overlapping_sketches->{$overlap}++;
					push(@{$hash->{$sketch}},[$id,2]);
				} else {
					$tooshort++;
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
print "Too short:".$tooshort."\n";
foreach my $key (keys(%{$overlapping_sketches})) {
	print $key."\t".$overlapping_sketches->{$key}."\n";
}
my $stats = {};
print "Scan start:".time()."\n";
for (my $i=0; $i < @{$genomelist}; $i++) { 
	print "Scanning:".$i." of ".$totalgenomes."\n";
	my $lines = &LOADFILE("/vol/patric3/downloads/genomes/".$genomelist->[$i]."/".$genomelist->[$i].".PATRIC.faa");
	my $id;
	my $func;
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/^>([^\s^\t]+)/) {
			my $newid = $1;
			my $hitcount = 0;
			my $one = 0;
			my $two = 0;
			my $zero = 0;
			my $hits = {};
			my $specieshash = {};
			my $functions = {};
			if (defined($id)) {
				if (50+$size < length($seq)) {
					for (my $k=0; $k < (length($seq)-$size-1); $k++) {
						my $query = substr($seq,$k,$size);
						if ($size > 30) {
							$query = Digest::MD5::md5_hex($query);
						}
						if (defined($hash->{$query})) {
							for (my $m=0; $m < @{$hash->{$query}}; $m++) {
								if (!defined($hits->{$hash->{$query}->[$m]->[0]})) {
									$hits->{$hash->{$query}->[$m]->[0]} = 0;
								}
								$hitcount++;
								$hits->{$hash->{$query}->[$m]->[0]}++;
								if (defined($funchash->{$hash->{$query}->[$m]->[0]})) {
									my $func = $funchash->{$hash->{$query}->[$m]->[0]};
									if (!defined($functions->{$func})) {
										$functions->{$func} = {
											"0" => 0,
											"1" => 0,
											"2" => 0,
											"g" => 0	
										};
									}
									if ($hash->{$query}->[$m]->[1] eq "0") {
										$zero++;
									}
									if ($hash->{$query}->[$m]->[1] eq "1") {
										$one++;
									}
									if ($hash->{$query}->[$m]->[1] eq "2") {
										$two++;
									}
									$functions->{$func}->{$hash->{$query}->[$m]->[1]}++;
									if ($hits->{$hash->{$query}->[$m]->[0]} == 1) {
										$functions->{$func}->{g}++;
									}
								}
								if ($hash->{$query}->[$m]->[0] =~ m/\|(\d+\.\d+)\./) {
									my $species = $genomes->{$1}->{n};
									my $array = [split(/\s/,$species)];
									if (!defined($specieshash->{$array->[0]})) {
										$specieshash->{$array->[0]} = {
											"0" => 0,
											"1" => 0,
											"2" => 0,
											"g" => 0
										};
									}
									$specieshash->{$array->[0]}->{$hash->{$query}->[$m]->[1]}++;
									if ($hits->{$hash->{$query}->[$m]->[0]} == 1) {
										$specieshash->{$array->[0]}->{g}++;
									}
								}	
							}
						}
					}
					my $genecount = keys(%{$hits});
					my $genuscount = keys(%{$specieshash});
					my $funccount = keys(%{$functions});
					if (!defined($stats->{middist}->{$one})) {
						$stats->{middist}->{$one} = 0;
					}
					$stats->{middist}->{$one}++;
					if (!defined($stats->{enddist}->{$two})) {
						$stats->{enddist}->{$two} = 0;
					}
					$stats->{enddist}->{$two}++;
					if (!defined($stats->{begdist}->{$zero})) {
						$stats->{begdist}->{$zero} = 0;
					}
					$stats->{begdist}->{$zero}++;
					if (!defined($stats->{sketchdist}->{$hitcount})) {
						$stats->{sketchdist}->{$hitcount} = 0;
					}
					$stats->{sketchdist}->{$hitcount}++;
					if (!defined($stats->{genedist}->{$genecount})) {
						$stats->{genedist}->{$genecount} = 0;
					}
					$stats->{genedist}->{$genecount}++;
					if (!defined($stats->{genusdist}->{$genuscount})) {
						$stats->{genusdist}->{$genuscount} = 0;
					}
					$stats->{genusdist}->{$genuscount}++;
					if (!defined($stats->{funcdist}->{$funccount})) {
						$stats->{funcdist}->{$funccount} = 0;
					}
					$stats->{funcdist}->{$funccount}++;
					my $list = [sort { $specieshash->{$a}->{g} <=> $specieshash->{$b}->{g} } keys(%{$specieshash})];
					for (my $m=0; $m < @{$list}; $m++) {
						if (defined($list->[$m])) {
							my $fract = $specieshash->{$list->[$m]}->{g}/$genecount;
							$fract = floor($fract*10);
							if (!defined($stats->{genusfract}->[$m]->[$fract])) {
								$stats->{genusfract}->[$m]->[$fract] = 0;
							}
							$stats->{genusfract}->[$m]->[$fract]++;
						}
					}
					$list = [sort { $functions->{$a} <=> $functions->{$b} } keys(%{$functions})];
					for (my $m=0; $m < @{$list}; $m++) {
						if (defined($list->[$m])) {
							my $fract = $functions->{$list->[$m]}->{g}/$genecount;
							$fract = floor($fract*10);
							if (!defined($stats->{funcfract}->[$m]->[$fract])) {
								$stats->{funcfract}->[$m]->[$fract] = 0;
							}
							$stats->{funcfract}->[$m]->[$fract]++;
						}
					}			
				}
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


my $output = &compute_stats($stats,$hash,$size,1);
&PRINTFILE($path."/stats.json",[&TOJSON($stats,1)]);

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