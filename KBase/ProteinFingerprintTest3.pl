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
	print "Sketching:".$i." of ".$totalgenomes."\t".time()."\n";
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
	#if ($count > 1000) {
	#	last;
	#}
}
print "Count:".$count."\n";
print "Too short:".$tooshort."\n";
foreach my $key (keys(%{$overlapping_sketches})) {
	print $key."\t".$overlapping_sketches->{$key}."\n";
}

my $filelist = [qw(
GO-all
GO-ERR2162200
GO-ERR2162201
GO-ERR2162202
GO-ERR2162203
GO-ERR2162204
GO-ERR2162205
GO-ERR2162206
GO-ERR2162207
GO-ERR2162208
GO-ERR2162209
GO-ERR2162210
GO-ERR2162211
GO-ERR2162212
GO-ERR2162213
GO-ERR2162214
GO-ERR2162215
GO-ERR2162216
GO-ERR2162217
GO-ERR2162218
GO-ERR2162219
GO-ERR2162220
GO-ERR2162221
GO-ERR2162222
GO-ERR2162223
GO-ERR2162224
)];
my $output = {};
for (my $i=0; $i < @{$filelist}; $i++) {
	print "Gene scan ".$filelist->[$i]."\n";
	my $features = &FROMJSON(join("\n",@{&LOADFILE($path."/all_genes/".$filelist->[$i].".json")}));
	foreach my $contigid (keys(%{$features})) {
		for (my $j=0; $j < @{$features->{$contigid}}; $j++) {
			my $seq = $features->{$contigid}->[$j]->{protein_translation};
			my $high = $features->{$contigid}->[$j]->{location}->[0]->[1] + $features->{$contigid}->[$j]->{location}->[0]->[3];
			my $low = $features->{$contigid}->[$j]->{location}->[0]->[1];
			if ($features->{$contigid}->[$j]->{location}->[0]->[2] eq "-") {
				$low = $features->{$contigid}->[$j]->{location}->[0]->[1] - $features->{$contigid}->[$j]->{location}->[0]->[3];
				$high = $features->{$contigid}->[$j]->{location}->[0]->[1];
			}
			if (!defined($output->{$contigid}->{lowest}) || $output->{$contigid}->{lowest} > $low) {
				$output->{$contigid}->{lowest} = $low;
			}
			if (!defined($output->{$contigid}->{highest}) || $output->{$contigid}->{highest} < $high) {
				$output->{$contigid}->{highest} = $high;
			}
			if ($size < length($seq)) {
				for (my $k=0; $k < (length($seq)-$size-1); $k++) {
					my $query = substr($seq,$k,$size);
					if ($size > 30) {
						$query = Digest::MD5::md5_hex($query);
					}
					my $pointer = $hash->{$query};
					if (defined($pointer)) {
						for (my $m=0; $m < @{$pointer}; $m++) {
							my $function = "unknown";
							if (defined($funchash->{$pointer->[$m]->[0]})) {
								$function = $funchash->{$pointer->[$m]->[0]};
							}
							my $species = "unknown";
							if ($pointer->[$m]->[0] =~ m/\|(\d+\.\d+)\./) {
								$species = $genomes->{$1}->{n};
							}
							$output->{$contigid}->{genes}->{$features->{$contigid}->[$j]->{id}}->{$pointer->[$m]->[0]} = {
								f => $function,
								s => $species
							};	
						}
					}
				}		
			}
		}
	}
}
for (my $i=0; $i < @{$filelist}; $i++) {
	print "Contig scan ".$filelist->[$i]."\n";
	my $lines = &LOADFILE($path."/all_genes/".$filelist->[$i].".json");
	my $id;
	my $func;
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/^>([^\s^\t]+)/) {
			my $newid = $1;
			if (defined($id)) {
				if (defined($output->{$id})) {
					if ($output->{$id}->{lowest} > 3*$size) {
						&CheckContigSequence(substr($seq,0,$output->{$id}->{lowest}),$output,"head",$id);
					}
					if (length($seq) > $output->{$id}->{highest}+3*$size) {
						&CheckContigSequence(substr($seq,$output->{$id}->{highest}),$output,"tail",$id);
					}
				} else {
					&CheckContigSequence($seq,$output,"entire",$id);
				}
				$id = $newid;
				$seq = "";
			} else {
				$seq .= $lines->[$j];
			}
		}
	}	
}

&PRINTFILE($path."/output-".$size.".json",[&TOJSON($output,1)]);

sub CheckContigSequence {
    my ($seq,$outhash,$type,$contigid) = @_;
	my $protseq = GUSTPlus::gustoenv::translate_sequence($seq,1);
	&ScanProteinForHits($protseq,$outhash,$type."_F1",$contigid);
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($seq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F2",$contigid);
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($seq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F3",$contigid);
	my $revseq = GUSTPlus::gustoenv::reverse_sequence($seq);
	$protseq = GUSTPlus::gustoenv::translate_sequence($revseq,1);
	&ScanProteinForHits($protseq,$outhash,$type."_R1",$contigid);
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($revseq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R2",$contigid);
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($revseq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R3",$contigid);	
}

sub ScanProteinForHits {
	my ($seq,$outhash,$type,$contigid,$id) = @_;
	if ($size < length($seq)) {
		for (my $k=0; $k < (length($seq)-$size-1); $k++) {
			my $query = substr($seq,$k,$size);
			if ($size > 30) {
				$query = Digest::MD5::md5_hex($query);
			}
			my $pointer = $hash->{$query};
			if (defined($pointer)) {
				for (my $m=0; $m < @{$pointer}; $m++) {
					my $function = "unknown";
					if (defined($funchash->{$pointer->[$m]->[0]})) {
						$function = $funchash->{$pointer->[$m]->[0]};
					}
					my $species = "unknown";
					if ($pointer->[$m]->[0] =~ m/\|(\d+\.\d+)\./) {
						$species = $genomes->{$1}->{n};
					}
					if ($type eq "protein") {
						$outhash->{$contigid}->{genes}->{$id}->{$pointer->[$m]->[0]} = {
							f => $function,
							s => $species
						};
					} else {
						$outhash->{$contigid}->{raw}->{$pointer->[$m]->[0]} = {
							type => $type,
							f => $function,
							s => $species
						};
					}	
				}
			}
		}		
	}
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