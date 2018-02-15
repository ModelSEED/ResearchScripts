#!/usr/bin/perl -w

use strict;
use JSON::XS;
use JSON;

local $| = 1;

my $path = "/disks/p3dev3/PATRICDump/";
my $kmer_dir = "/vol/patric3/fams/2017-0701/kmers";
my $url = "http://spruce:6100";
my $core_dir = "/disks/p3dev3/code/Genome_quality/Genus_Core_PGFs";
my $fasta_dir = "/vol/patric3/downloads/genomes/";
my $output_dir = "/disks/p3dev3/PATRICDump/";

my $genomes = &FROMJSON(join("\n",@{&LOADFILE($path."/AllPatricGenomes.txt")}));

my $success = 0;
my $fail = 0;
foreach my $genome (keys(%{$genomes})) {
	my $genus;
	if (-e $fasta_dir.$genome."/".$genome.".PATRIC.faa") {
		if ($genomes->{$genome}->{n} =~ m/^(.+)\s/) {
			$genus = $1;
		}
		
		system("perl /disks/p3dev3/code/Genome_quality/core_PGF_genome_quality.pl -c ".$core_dir." -d ".$kmer_dir." -u ".$url." -g Escherichia -p ".$fasta_dir.$genome."/".$genome.".PATRIC.faa > ".$output_dir."compout");
		
		if (-e $output_dir."compout") {
			my $output = &LOADFILE($output_dir."compout");
			if (defined($output->[1])) {
				print $genome."\t".$output->[1]."\n";
				my $array = [split(/\t/,$output->[1])];
				if (defined($array->[7])) {
					$genomes->{$genome}->{ma} = $array->[0];
					$genomes->{$genome}->{mf} = $array->[1];
					$genomes->{$genome}->{fcom} = $array->[2];
					$genomes->{$genome}->{fdup} = $array->[3];
					$genomes->{$genome}->{fmiss} = $array->[4];
					$genomes->{$genome}->{fshort} = $array->[5];
					$genomes->{$genome}->{flong} = $array->[6];
					$genomes->{$genome}->{ratio} = $array->[7];
					$success++;
				}
			}
			unlink($output_dir."compout");
		} else {
			print $genome." failed!\n";
			$fail++;
		}
	} else {
		print $genome." failed!\n";
		$fail++;
	}
}

&PRINTFILE($path."/AllPatricGenomesNew.txt",[&TOJSON($genomes,1)]);

print "\n\nSuccess:".$success."\nFail:".$fail."\n";

sub PRINTFILE {
    my ($filename,$arrayRef) = @_;
    open ( my $fh, ">", $filename);
    foreach my $Item (@{$arrayRef}) {
    	print $fh $Item."\n";
    }
    close($fh);
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