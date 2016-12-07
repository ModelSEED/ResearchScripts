#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use JSON;
use Bio::KBase::workspaceService::Helpers qw(auth get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::Helpers qw(get_fba_client runFBACommand universalFBAScriptCode );

my $mediaCount = 0;
my $mediaHash;
my $mediaList;
my $rxnData;
my $essMediaPattern;
my $index = 3;
for (my $i=0;$i < 50; $i++) {
	my $index = 3+$i;
	my $output = runFBACommand({
		fbas => ["kb|g.20848.fbamdl.17.fba.".$index],
		workspaces => ["phenotypeDemo"],
		auth => auth()
	},"get_fbas",{showerror => 1});
	if (defined($output->[0]->{formulation}->{media})) {
		my $media = $output->[0]->{formulation}->{media};
		if (defined($output->[0]->{reactionFluxes})) {
			for (my $i=0;$i < @{$output->[0]->{reactionFluxes}}; $i++) {
				my $rxn = $output->[0]->{reactionFluxes}->[$i];
				my $maxflux = $rxn->[4];
				my $minflux = $rxn->[5];
				my $id = $rxn->[0];
				if ($minflux > 0) {
					$rxnData->{$id}->{$media} = 1;
				} elsif ($maxflux < 0) {
					$rxnData->{$id}->{$media} = 1;
				}
			}
		}
		if (!defined($mediaHash->{$media})) {
			$mediaHash->{$media} = 1;
			$mediaCount++;
		}
	}
}

foreach my $rxn (keys(%{$rxnData})) {
	my $mediaList = join(";",sort(keys(%{$rxnData->{$rxn}})));
	$essMediaPattern->{$mediaList}->{reactions}->{$rxn} = 1;
}
my $serv = get_ws_client();
my $output = $serv->get_object({
	id => "kb|g.20848.phenos.21",
	type => "PhenotypeSet",
	workspace => "phenotypeDemo",
	auth => auth()
});
my $data = $output->{data};
my $geneData;
for (my $i=0;$i < @{$data->{phenotypes}}; $i++) {
	my $phenotype = $data->{phenotypes}->[$i];
	if (defined($phenotype->[0]->[0]) && $phenotype->[4] eq "0") {
		$geneData->{$phenotype->[0]->[0]}->{$phenotype->[1]} = 1;
	}
}
foreach my $gene (keys(%{$geneData})) {
	my $mediaList = join(";",sort(keys(%{$geneData->{$gene}})));
	$essMediaPattern->{$mediaList}->{genes}->{$gene} = 1;
}		

my $model = runFBACommand({
	models => ["kb|g.20848.fbamdl.17"],
	workspaces => ["phenotypeDemo"],
	auth => auth()
},"get_models",{showerror => 1});
my $rxnGenes;
my $geneRxns;
if (defined($model->[0]->{reactions})) {
	for (my $i=0;$i < @{$model->[0]->{reactions}}; $i++) {
		my $rxn = $model->[0]->{reactions}->[$i];
		if (defined($rxn->{features})) {
			for (my $j=0;$j < @{$rxn->{features}}; $j++) {
				$rxnGenes->{$rxn->{id}}->{$rxn->{features}->[$j]} = 1;
				$geneRxns->{$rxn->{features}->[$j]}->{$rxn->{id}} = 1;
			}
		}
	}
}
$serv = get_ws_client();
$output = $serv->get_object({
	id => "kb|g.20848.phenos.21.simulation",
	type => "PhenotypeSimulationSet",
	workspace => "phenotypeDemo",
	instance => 2,
	auth => auth()
});
print "Essential media list\tUnmodeled genes\tModeled genes\tReactions with genes\tReactions without genes\n";
foreach my $mediaList (keys(%{$essMediaPattern})) {
	my $modGenes = [];
	my $unmodGenes = [];
	foreach my $gene (keys(%{$essMediaPattern->{$mediaList}->{genes}})) {
		if (defined($geneRxns->{$gene})) {
			push(@{$modGenes},$gene);
		} else {
			push(@{$unmodGenes},$gene);
		}
	}
	my $geneRxn = [];
	my $noGeneRxn = [];
	foreach my $rxn (keys(%{$essMediaPattern->{$mediaList}->{reactions}})) {
		if (defined($rxnGenes->{$rxn})) {
			push(@{$geneRxn},$rxn);
		} else {
			push(@{$noGeneRxn},$rxn);
		}
	}
	print $mediaList."\t".join(";",@{$unmodGenes})."\t".join(";",@{$modGenes})."\t".join(";",@{$geneRxn})."\t".join(";",@{$noGeneRxn})."\n";		
}
print "\n";
my $missingIsozymes;
my $badComplex;
my $extraIsozyme;
my $essRxnEssGenes;
$data = $output->{data};
for (my $i=0;$i < @{$data->{phenotypeSimulations}}; $i++) {
	my $simulation = $data->{phenotypeSimulations}->[$i];
	my $phenotype = $simulation->[0];
	if (defined($phenotype->[0]->[0])) {
		my $media = $phenotype->[1];
		my $geneko = $phenotype->[0]->[0];
		my $growthFraction = $simulation->[2];
		my $class = $simulation->[3];
		if ($class eq "FP") {
			#Scan through reactions associated with gene
			foreach my $rxn (keys(%{$geneRxns->{$geneko}})) {
				#Check if the reaction is essential in the specified media
				if (defined($rxnData->{$rxn}->{$media})) {
					#Check if the reaction has multiple genes
					if (keys(%{$rxnGenes->{$rxn}}) > 1) {
						#Reaction appears to have extra isozymes
						$extraIsozyme->{$rxn}->{$geneko}->{$media} = 1;
					} else {
						#Reaction is essential and mapped to essential gene, indicating a problem
						$essRxnEssGenes->{$rxn}->{$geneko}->{$media} = 1;
					}
				}
			}
		} elsif ($class eq "FN") {
			#Scan through reactions associated with gene
			foreach my $rxn (keys(%{$geneRxns->{$geneko}})) {
				#Check if the reaction is essential in the specified media
				if (defined($rxnData->{$rxn}->{$media})) {
					#Check if the reaction has multiple genes
					if (keys(%{$rxnGenes->{$rxn}}) > 1) {
						#Reaction appears to have extra isozymes
						$badComplex->{$rxn}->{$geneko}->{$media} = 1;
					} else {
						#Reaction is essential and mapped to essential gene, indicating missing isozyme
						$missingIsozymes->{$rxn}->{$geneko}->{$media} = 1;
					}
				}
			}
		}
	}
}
print "Reaction with missing isozyme\tCurrent gene\tProblematic media\n";
foreach my $rxn (keys(%{$missingIsozymes})) {
	foreach my $gene (keys(%{$missingIsozymes->{$rxn}})) {
		print $rxn."\t".$gene."\t".join(";",keys(%{$missingIsozymes->{$rxn}->{$gene}}))."\n";
	}
}
print "\n";
print "Essential rxn\tEssential gene\tProblematic media\n";
foreach my $rxn (keys(%{$essRxnEssGenes})) {
	foreach my $gene (keys(%{$essRxnEssGenes->{$rxn}})) {
		print $rxn."\t".$gene."\t".join(";",keys(%{$essRxnEssGenes->{$rxn}->{$gene}}))."\n";
	}
}
print "\n";
print "Essential rxn\tEssential genes\tNonessential genes\n";
foreach my $rxn (keys(%{$badComplex})) {
	my $essentialGenes = [];
	my $nonEssentialGenes = [];
	foreach my $gene (keys(%{$badComplex->{$rxn}})) {
		push(@{$essentialGenes},$gene."(".join(";",keys(%{$badComplex->{$rxn}->{$gene}})).")");
	}
	foreach my $gene (keys(%{$rxnGenes->{$rxn}})) {
		if (!defined($badComplex->{$rxn}->{$gene})) {
			push(@{$nonEssentialGenes},$gene);
		}
	}
	print $rxn."\t".join(";",@{$essentialGenes})."\t".join(";",@{$nonEssentialGenes})."\n";
}
print "\n";
print "Essential rxn\tEssential genes\tNonessential genes\n";
foreach my $rxn (keys(%{$extraIsozyme})) {
	my $essentialGenes = [];
	my $nonEssentialGenes = [];
	foreach my $gene (keys(%{$extraIsozyme->{$rxn}})) {
		push(@{$essentialGenes},$gene."(".join(";",keys(%{$extraIsozyme->{$rxn}->{$gene}})).")");
	}
	foreach my $gene (keys(%{$rxnGenes->{$rxn}})) {
		if (!defined($extraIsozyme->{$rxn}->{$gene})) {
			push(@{$nonEssentialGenes},$gene);
		}
	}
	print $rxn."\t".join(";",@{$essentialGenes})."\t".join(";",@{$nonEssentialGenes})."\n";
}
print "\n";