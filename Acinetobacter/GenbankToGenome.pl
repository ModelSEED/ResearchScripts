#!/usr/bin/perl

use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

$|=1;

my $filename = $ARGV[0];

my $linearray = [];
open(my $fh, "<", $filename);
while (my $line = <$fh>) {
	chomp($line);
	push(@{$linearray},$line);
}
close($fh);

my $genome = {
	scientific_name => "Acinetobacter baumannii 5075",
	features => [],
	domain => "Bacteria",
	genetic_code => 11,
	id => "Acinetobacter_baumannii_5075",
	contigset_ref => "8834/8/2",
	gc_content => 0.5,
	dna_size => 4051011,
	num_contigs => 1,
	contig_lengths => [4051011],
	contig_ids => ["5075"],
	source => "http://www.gs.washington.edu/labs/manoil/baumannii.htm",
	source_id => "Acinetobacter_baumannii_5075",
	taxonomy => "Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter baumannii AYE",
	publications => [],
	close_genomes => [],
	analysis_events => []
};
my $hash;
for (my $i=0; $i < @{$linearray}; $i++) {
	if ($linearray->[$i] =~ m/^\s\s\s\s\s([A-Za-z]+)\s+(\d+)\.\.(\d+)/ || $linearray->[$i] =~ m/^\s\s\s\s\s([A-Za-z]+)\s+complement\((\d+)\.\.(\d+)\)/) {
		my $start = $2;
		my $stop = $3;
		my $gene = {
			type => $1,
			id => undef,
			location => [
				[
					"5075",
					$start,
					"+",
					($stop-$start)
				]
			]
		};
		if ($linearray->[$i] =~ m/complement/) {
			$gene->{location}->[0]->[2] = "-";
			$gene->{location}->[0]->[1] = $stop+0;
		}
		$gene->{location}->[0]->[3] += 0;
		$gene->{location}->[0]->[1] += 0;
		$i++;
		while ($linearray->[$i] !~ /^\s\s\s\s\s[A-Za-z]/ && $i < @{$linearray}) {
			if ($linearray->[$i] =~ m/locus_tag=\"(.+)\"$/) {
				$gene->{id} = $1;
			} elsif ($linearray->[$i] =~ m/product=\"(.+)\"$/) {
				$gene->{function} = $1;
			} elsif ($linearray->[$i] =~ m/translation=\"([^\"]+)/) {
				$gene->{protein_translation} = $1;
				while($linearray->[$i+1] =~ m/^\s+([A-Z]+)\"*$/) {
					$gene->{protein_translation} .= $1;
					$i++;
				}
			}
			$i++;
		}
		$i--;
		if (defined($gene->{id}) && defined($gene->{function})) {
			$hash->{$gene->{id}} = $gene;
		}
	}
}
foreach my $key (keys(%{$hash})) {
	#print $key."\t".$hash->{$key}->{type}."\t".join("\t",@{$hash->{$key}->{location}->[0]})."\t".$hash->{$key}->{function}."\t".$hash->{$key}->{protein_translation}."\n";
	push(@{$genome->{features}},$hash->{$key});
}

save_workspace_object("chenry:1436155557602/Acinetobacter_baumannii_5075",$genome,"KBaseGenomes.Genome");

1;