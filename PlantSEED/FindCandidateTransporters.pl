#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;
use Data::Dumper;

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

print "Loading genome list!\n";
my $genomelist;
open(my $fhh, "<", $directory."GenomeList.txt");
while (my $line = <$fhh>) {
	chomp($line);
	push(@{$genomelist},$line);
}
close($fhh);

print "Loading blast data!\n";
open(my $fhhh, "<", $directory."BlastOutput.txt");
my $genehash;
while (my $line = <$fhhh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$genehash->{$array->[0]} = [split(/;/,$array->[1])];
}
close($fhhh);

my $alreadydone = {};
my $candidategenes;
print "Searching for transporter candidates!\n";
for (my $i=0; $i < @{$genomelist}; $i++) {
#for (my $i=0; $i < 3; $i++) {
	my $genome = $genomelist->[$i];
	if (!defined($alreadydone->{$genome})) {
		print "Loading genome ".$genome."\n";
		my $genomeHash = $sapsvr->all_features({
			-ids => [$genome],
		});
		print Data::Dumper->Dump([$genomeHash]);
		my $lochash = $sapsvr->fid_locations({
			-ids => $genomeHash->{$genome}
		});
		my $functions = $sapsvr->ids_to_functions({-ids => $genomeHash->{$genome}});
		for (my $j=0; $j < @{$genomeHash->{$genome}}; $j++) {
			#print $genomeHash->{$genome}->[$j]."\t".$lochash->{$genomeHash->{$genome}->[$j]}->[0]."\t".$functions->{$genomeHash->{$genome}->[$j]}."\n";
		}
		my $contigs = {};
		my $geneindex = {};
		foreach my $gene (keys(%{$lochash})) {
			if ($lochash->{$gene}->[0] =~ m/(.+)_(\d+)_(\d+)$/) {
				push(@{$contigs->{$1}},$gene);
				$geneindex->{$gene} = $2;
			} elsif ($lochash->{$gene}->[0] =~ m/(.+)_(\d+)[\+-](\d+)$/) {
				push(@{$contigs->{$1}},$gene);
				$geneindex->{$gene} = $2;
			}
		}
		foreach my $gene (keys(%{$functions})) {
			$functions->{$gene} =~ s/\s*#.+//;
			$functions->{$gene} = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$functions->{$gene}))];
			for (my $j=0; $j < @{$functions->{$gene}}; $j++) {
				$functions->{$gene}->[$j] = lc($functions->{$gene}->[$j]);
				$functions->{$gene}->[$j] =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
				$functions->{$gene}->[$j] =~ s/\s//g;
				$functions->{$gene}->[$j] =~ s/\#.*$//g;
				$functions->{$gene}->[$j] =~ s/\(ec\)//g;
			}
		}
		foreach my $contig (keys(%{$contigs})) {
			$contigs->{$contig} = [sort { $geneindex->{$a} <=> $geneindex->{$b} } @{$contigs->{$contig}}];
			for (my $j=0;$j < @{$contigs->{$contig}}; $j++) {
				for (my $n=0; $n < @{$functions->{$contigs->{$contig}->[$j]}}; $n++) {
					if (defined($rolehash->{$functions->{$contigs->{$contig}->[$j]}->[$n]})) {
						my $start = $j-5;
						my $stop = $j+6;
						if ($start < 0) {
							$start = 0;
						}
						if ($stop > @{$contigs->{$contig}}) {
							$stop = @{$contigs->{$contig}};
						}
						for (my $k=$start;$k < $stop; $k++) {
							for (my $m=0;$m < @{$functions->{$contigs->{$contig}->[$k]}}; $m++) {
								if (!defined($rolehash->{$functions->{$contigs->{$contig}->[$k]}->[$m]})) {
									if (defined($genehash->{$contigs->{$contig}->[$k]})) {
										for (my $p=0; $p < @{$genehash->{$contigs->{$contig}->[$k]}}; $p++) {
											$candidategenes->{$genehash->{$contigs->{$contig}->[$k]}->[$p]}->{$contigs->{$contig}->[$k]} = 1;											
										}
									}
								}
							}
						}	
					}
				}
			}
		}
	}
}
print "Printing candidates:\n";
foreach my $candidate (keys(%{$candidategenes})) {
	print $candidate."\t".join(";",keys(%{$candidategenes->{$candidate}}))."\n";
}