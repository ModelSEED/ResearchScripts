#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];
my $window = 5;

my $sapsvr = ModelSEED::Client::SAP->new();

print "Loading vitamins!\n";
open(my $fh, "<", $directory."BVitaminRoles.txt");
my $rolehash = {};
while (my $line = <$fh>) {
	chomp($line);
	my $rolearray = [split(/\t/,$line)];
	if (defined($rolearray->[3])) {
		my $itemarray = [split(/;/,$rolearray->[3])];
		if (defined($itemarray->[2])) {
			$rolehash->{$itemarray->[2]} = $itemarray->[1];
			$rolehash->{$itemarray->[1]} = $itemarray->[1];
		}
	}
}
close($fh);

print "Loading blast file list!\n";
my $array;
open(my $fhh, "<", $directory."BlastFileList.txt");
while (my $line = <$fhh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fhh);

print "Loading blast data!\n";
my $genomehash;
my $plantgenehash;
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i." ".$array->[$i]."\n";
	open(my $fhhh, "<", "/homes/seaver/Projects/Plants_PubSEED_Sims/Blast_Results/".$array->[$i]);
	while (my $line = <$fhhh>) {
		chomp($line);
		my $blastdata = [split(/\t/,$line)];
		if ($blastdata->[1] =~ m/(fig\|\d+\.\d+)\./) {
			$genomehash->{$1} = 1;
			$plantgenehash->{$blastdata->[1]}->{$blastdata->[0]} = 1;
		}
	}
	close($fhhh);
}

my $alreadydone = {};
print "Loading genome data of ".keys(%{$genomehash})."!\n";
foreach my $genome (keys(%{$genomehash})) {
	if (!defined($alreadydone->{$genome})) {
		print "Loading genome ".$genome."\n";
		my $colocalizedgenes;
		my $genomeHash = $sapsvr->all_features({
			-ids => [$genome],
		});
		my $lochash = $sapsvr->fid_locations({
			-ids => $genomeHash->{$genome}
		});
		my $functions = $sapsvr->ids_to_functions({-ids => $genomeHash->{$genome}});
		foreach my $gene (keys(%{$functions})) {
			$functions->{$gene} = lc($functions->{$gene});
			$functions->{$gene} =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
			$functions->{$gene} =~ s/\s//g;
			$functions->{$gene} =~ s/\#.*$//g;
			$functions->{$gene} =~ s/\(ec\)//g;
		}
		my $contigs = {};
		my $geneindex = {};
		foreach my $gene (keys(%{$lochash})) {
			if ($lochash->{$gene} =~ m/(.+)_(\d+)_(\d+)$/) {
				push(@{$contigs->{$1}},$gene);
				$geneindex->{$gene} = $2;
			} elsif ($lochash->{$gene} =~ m/(.+)_(\d+)[\+-](\d+)$/) {
				push(@{$contigs->{$1}},$gene);
				$geneindex->{$gene} = $2;
			}
		}
		foreach my $contig (keys(%{$contigs})) {
			$contigs->{$contig} = [sort { $geneindex->{$a} <=> $geneindex->{$b} } @{$contigs->{$contig}}];
			for (my $i=0;$i < @{$contigs->{$contig}}; $i++) {
				if (defined($rolehash->{$functions->{$contigs->{$contig}->[$i]}})) {
					my $start = $i-5;
					my $stop = $i+6;
					if ($start < 0) {
						$start = 0;
					}
					if ($stop > @{$contigs->{$contig}}) {
						$stop = @{$contigs->{$contig}};
					}
					for (my $j=$start;$j < $stop; $j++) {
						if (!defined($rolehash->{$functions->{$contigs->{$contig}->[$j]}})) {
							if (!defined($plantgenehash->{$contigs->{$contig}->[$j]})) {
								$colocalizedgenes->{$contigs->{$contig}->[$j]} = 1;
							}
						}
					}
				}
			}
		}
		foreach my $gene (keys(%{$colocalizedgenes})) {
			print "Colocalized genes:".join(";",keys(%{$colocalizedgenes}))."\n";
		}
	}
}