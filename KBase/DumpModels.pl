#!/usr/bin/perl -w

use strict;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
$|=1;

my $directory = $ARGV[0];

open ( my $fh, ">", $directory."/ModelData.txt");
open ( my $fhh, ">", $directory."/ModelTable.txt");

print $fh "ID\tName\tSource ID\tGenome ID\tGenome source ID\tGenome name\tTaxonomy\tDomain\tSize\tContigs\n";
print $fhh "ID\tReaction ID\tDirections\tGenes\n";
my $list = get_ws_objects_list("KBasePublicModelsV4","KBaseFBA.FBAModel");
for (my $i=0; $i < @{$list}; $i++) {
	if ($list->[$i]->[1] =~ m/\.gf$/) {
		print "Getting model ".$list->[$i]->[1]."\n";
		(my $obj,my $meta) = get_workspace_object("KBasePublicModelsV4/".$list->[$i]->[1]);
		print "Getting genome ".$obj->{genome_ref}."\n";
		(my $genome,my $meta) = get_workspace_object($obj->{genome_ref});
		print "Printing data...\n";
		print $fh $obj->{id}."\t".$obj->{name}."\t".$obj->{source_id}."\t".$genome->{id}."\t".$genome->{source_id}."\t".$genome->{scientific_name}."\t".$genome->{taxonomy}."\t".$genome->{domain}."\t".$genome->{dna_size}."\t".$genome->{num_contigs}."\n";
		my $rxns = $obj->{modelreactions};
		for (my $j=0; $j < @{$rxns}; $j++) {
			my $prots = $rxns->[$j]->{modelReactionProteins};
			for (my $k=0; $k < @{$prots}; $k++) {
				my $subunits = $prots->[$k]->{modelReactionProteinSubunits};
				for (my $m=0; $m < @{$subunits}; $m++) {
					my $features = $subunits->[$m]->{feature_refs};
					$subunits->[$m] = join(" or" ,@{$features});
					if (@{$features} > 1) {
						$subunits->[$m] = "(".$subunits->[$m].")";
					} elsif (@{$features} == 0) {
						$subunits->[$m] = "Unknown";
					}
				}
				$prots->[$k] = join(" and ",@{$subunits});
				if (@{$subunits} > 1) {
					$prots->[$k] = "(".$prots->[$k].")";
				}
			}
			my $genes = join(" or ",@{$prots});
			if (@{$prots} > 1) {
				$genes = "(".$genes.")";
			}
			if (length($genes) == 0) {
				$genes = "Unknown";
			}
			print $fhh $obj->{id}."\t".$rxns->[$j]->{id}."\t".$rxns->[$j]->{direction}."\t".$genes."\n";
		}
	}
}

close($fh);
close($fhh);