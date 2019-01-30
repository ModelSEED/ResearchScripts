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
		$output->{$filelist->[$i]}->{$contigid} = {
			genes => {},
			raw => {},
			genehits => 0,
			rawhits => 0
		};
		for (my $j=0; $j < @{$features->{$contigid}}; $j++) {
			$output->{$filelist->[$i]}->{$contigid}->{genes}->{$features->{$contigid}->[$j]->{id}} = {
				s => {},
				f => {},
				count => 0
			};
			my $seq = $features->{$contigid}->[$j]->{protein_translation};
			my $high = $features->{$contigid}->[$j]->{location}->[0]->[1] + $features->{$contigid}->[$j]->{location}->[0]->[3];
			my $low = $features->{$contigid}->[$j]->{location}->[0]->[1];
			if ($features->{$contigid}->[$j]->{location}->[0]->[2] eq "-") {
				$low = $features->{$contigid}->[$j]->{location}->[0]->[1] - $features->{$contigid}->[$j]->{location}->[0]->[3];
				$high = $features->{$contigid}->[$j]->{location}->[0]->[1];
			}
			if (!defined($output->{$filelist->[$i]}->{$contigid}->{lowest}) || $output->{$filelist->[$i]}->{$contigid}->{lowest} > $low) {
				$output->{$filelist->[$i]}->{$contigid}->{lowest} = $low;
			}
			if (!defined($output->{$filelist->[$i]}->{$contigid}->{highest}) || $output->{$filelist->[$i]}->{$contigid}->{highest} < $high) {
				$output->{$filelist->[$i]}->{$contigid}->{highest} = $high;
			}
			&ScanProteinForHits($seq,$output->{$filelist->[$i]},"protein",$contigid,$features->{$contigid}->[$j]->{id});
			if ($output->{$filelist->[$i]}->{$contigid}->{genes}->{$features->{$contigid}->[$j]->{id}}->{count} > 0) {
				$output->{$filelist->[$i]}->{$contigid}->{genehits}++;
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
				if (defined($output->{$filelist->[$i]}->{$id})) {
					if ($output->{$filelist->[$i]}->{$id}->{lowest} > 3*$size) {
						$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence(substr($seq,0,$output->{$filelist->[$i]}->{$id}->{lowest}),$output->{$filelist->[$i]},"head",$id);
					}
					if (length($seq) > $output->{$filelist->[$i]}->{$id}->{highest}+3*$size) {
						$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence(substr($seq,$output->{$filelist->[$i]}->{$id}->{highest}),$output->{$filelist->[$i]},"tail",$id);
					}
				} else {
					$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence($seq,$output->{$filelist->[$i]},"entire",$id);
				}
				$id = $newid;
				$seq = "";
			} else {
				$seq .= $lines->[$j];
			}
		}
	}	
}

foreach my $id (keys(%{$output})) {
	my $genespeclines = ["Contig id\tGene id\tSpecies\tCount\tFraction"];
	my $genefunclines = ["Contig id\tGene id\tFunction\tCount\tFraction"];
	my $cmlines = ["Contig id\tGenes\tGeneHits\tRawHits\tHigh\tLow"];
	my $rawspeclines = ["Contig id\tType\tSpecies\tCount\tFraction"];
	my $rawfunclines = ["Contig id\tType\tFunction\tCount\tFraction"];
	foreach my $contigid (keys(%{$output->{$id}})) {
		my $genecount = keys(%{$output->{$id}->{$contigid}});
		my $high = "-";
		my $low = "-";
		if (defined($output->{$id}->{$contigid}->{lowest})) {
			$low = $output->{$id}->{$contigid}->{lowest};
		}
		if (defined($output->{$id}->{$contigid}->{lowest})) {
			$high = $output->{$id}->{$contigid}->{highest};
		}
		push(@{$cmlines},$contigid."\t".$genecount."\t".$output->{$id}->{$contigid}->{genehits}."\t".$output->{$id}->{$contigid}->{rawhits}."\t".$output->{$id}->{$contigid}->{rawhits}."\t".$high."\t".$low);
		foreach my $gene (keys(%{$output->{$id}->{$contigid}->{genes}})) {
			foreach my $species (keys(%{$output->{$id}->{$contigid}->{genes}->{$gene}->{s}})) {
				my $count = $output->{$id}->{$contigid}->{genes}->{$gene}->{s}->{$species};
				my $fraction = $count/$output->{$id}->{$contigid}->{genes}->{$gene}->{count};
				push(@{$genespeclines},$contigid."\t".$gene."\t".$species."\t".$count."\t".$fraction);
			}
			foreach my $function (keys(%{$output->{$id}->{$contigid}->{genes}->{$gene}->{f}})) {
				my $count = $output->{$id}->{$contigid}->{genes}->{$gene}->{f}->{$function};
				my $fraction = $count/$output->{$id}->{$contigid}->{genes}->{$gene}->{count};
				push(@{$genefunclines},$contigid."\t".$gene."\t".$function."\t".$count."\t".$fraction);
			}
		}
		foreach my $gene (keys(%{$output->{$id}->{$contigid}->{raw}})) {
			foreach my $species (keys(%{$output->{$id}->{$contigid}->{raw}->{$gene}->{s}})) {
				my $count = $output->{$id}->{$contigid}->{raw}->{$gene}->{s}->{$species};
				my $fraction = $count/$output->{$id}->{$contigid}->{raw}->{$gene}->{count};
				push(@{$rawspeclines},$contigid."\t".$gene."\t".$species."\t".$count."\t".$fraction);
			}
			foreach my $function (keys(%{$output->{$id}->{$contigid}->{raw}->{$gene}->{f}})) {
				my $count = $output->{$id}->{$contigid}->{raw}->{$gene}->{f}->{$function};
				my $fraction = $count/$output->{$id}->{$contigid}->{raw}->{$gene}->{count};
				push(@{$rawfunclines},$contigid."\t".$gene."\t".$function."\t".$count."\t".$fraction);
			}
		}
	}
	&PRINTFILE($path."/".$size."_".$id.".contig_meta",$cmlines);
	&PRINTFILE($path."/".$size."_".$id.".gene_spec",$genespeclines);
	&PRINTFILE($path."/".$size."_".$id.".gene_func",$genefunclines);
	&PRINTFILE($path."/".$size."_".$id.".contig_spec",$rawspeclines);
	&PRINTFILE($path."/".$size."_".$id.".contig_func",$rawfunclines);
}

sub CheckContigSequence {
    my ($seq,$outhash,$type,$contigid) = @_;
	my $count = 0;
	my $protseq = GUSTPlus::gustoenv::translate_sequence($seq,1);
	&ScanProteinForHits($protseq,$outhash,$type."_F1",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F1"}) && $outhash->{$contigid}->{raw}->{$type."_F1"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F1"}->{count};
	}
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($seq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F2",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F2"}) && $outhash->{$contigid}->{raw}->{$type."_F2"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F2"}->{count};
	}
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($seq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F3",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F3"}) && $outhash->{$contigid}->{raw}->{$type."_F3"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F3"}->{count};
	}
	my $revseq = GUSTPlus::gustoenv::reverse_sequence($seq);
	$protseq = GUSTPlus::gustoenv::translate_sequence($revseq,1);
	&ScanProteinForHits($protseq,$outhash,$type."_R1",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_R1"}) && $outhash->{$contigid}->{raw}->{$type."_R1"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R1"}->{count};
	}
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($revseq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R2",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_R2"}) && $outhash->{$contigid}->{raw}->{$type."_R2"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R2"}->{count};
	}
	$protseq = GUSTPlus::gustoenv::translate_sequence(substr($revseq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R3",$contigid);	
	if (defined($outhash->{$contigid}->{raw}->{$type."_R3"}) && $outhash->{$contigid}->{raw}->{$type."_R3"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R3"}->{count};
	}
	return $count;
}

sub ScanProteinForHits {
	my ($seq,$outhash,$type,$contigid,$id) = @_;
	if ($size < length($seq)) {
		my $genehash = {};
		for (my $k=0; $k < (length($seq)-$size-1); $k++) {
			my $query = substr($seq,$k,$size);
			if ($size > 30) {
				$query = Digest::MD5::md5_hex($query);
			}
			my $pointer = $hash->{$query};
			if (defined($pointer)) {
				for (my $m=0; $m < @{$pointer}; $m++) {
					my $new = 0;
					if (!defined($genehash->{$pointer->[$m]->[0]})) {
						$genehash->{$pointer->[$m]->[0]} = 0;
						$new = 1;
					}
					$genehash->{$pointer->[$m]->[0]}++;
					if ($new == 1) {
						my $function = "unknown";
						if (defined($funchash->{$pointer->[$m]->[0]})) {
							$function = $funchash->{$pointer->[$m]->[0]};
						}
						my $species = "unknown";
						if ($pointer->[$m]->[0] =~ m/\|(\d+\.\d+)\./) {
							$species = $genomes->{$1}->{n};
							my $array = [split(/\s/,$species)];
							if ($array->[1] eq "sp.") {
								$species = $array->[0]." ".$array->[2];
							} else {
								$species = $array->[0]." ".$array->[1];
							}
						}
						my $hash;
						if ($type eq "protein") {
							if (!defined($outhash->{$contigid}->{genes}->{$id})) {
								$outhash->{$contigid}->{genes}->{$id} = {s => {},f => {},count => 0};
							}
							$hash = $outhash->{$contigid}->{genes}->{$id};
						} else {
							if (!defined($outhash->{$contigid}->{raw}->{$type})) {
								$outhash->{$contigid}->{raw}->{$type} = {s => {},f => {},count => 0};
							}
							$hash = $outhash->{$contigid}->{raw}->{$type};
						}
						if (!defined($hash->{s}->{$species})) {
							$hash->{s}->{$species} = 0;
						}
						$hash->{s}->{$species}++;
						if (!defined($hash->{f}->{$function})) {
							$hash->{f}->{$function} = 0;
						}
						$hash->{f}->{$function}++;
						$hash->{count}++;
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