#!/usr/bin/perl -w

use strict;
use fba_tools::fba_toolsImpl;
$|=1;

Bio::KBase::kbaseenv::create_context_from_client_config();
my $impl = fba_tools::fba_toolsImpl->new();
my $fba = $impl->util_get_object("chenry:narrative_1526623886203/Maize_FBA_Complete");
my $rxns = $fba->FBAReactionBounds();
my $rxnhash;
for (my $i=0; $i < @{$rxns}; $i++) {
	$rxnhash->{$rxns->[$i]->modelreaction()->id()}->{complete} = {
		value => $rxns->[$i]->value(),
		class => $rxns->[$i]->class(),
	}
}
$fba = $impl->util_get_object("chenry:narrative_1526623886203/Maize_FBA_Heterotrophic");
$rxns = $fba->FBAReactionBounds();
for (my $i=0; $i < @{$rxns}; $i++) {
	$rxnhash->{$rxns->[$i]->modelreaction()->id()}->{heterotrophic} = {
		value => $rxns->[$i]->value(),
		class => $rxns->[$i]->class(),
	}
}
$fba = $impl->util_get_object("chenry:narrative_1526623886203/Maize_FBA_Autotophic");
$rxns = $fba->FBAReactionBounds();
for (my $i=0; $i < @{$rxns}; $i++) {
	$rxnhash->{$rxns->[$i]->modelreaction()->id()}->{autotrophic} = {
		value => $rxns->[$i]->value(),
		class => $rxns->[$i]->class(),
	}
}

my $keggtbl = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/code/fba_tools/data/KEGG_pathways");
my $kegghash = {};
for (my $i=1; $i < @{$keggtbl}; $i++) {
	my $items = [split(/\t/,$keggtbl->[$i])];
	$kegghash->{$items->[1]} = $items->[2];
}

my $rxndata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/code/fba_tools/data/Reactions.json")}));
my $rxndatahash = {};
for (my $i=0; $i < @{$rxndata}; $i++) {
	$rxndatahash->{$rxndata->[$i]->{id}} = $rxndata->[$i];
}

my $rxnec = {};
my $maizereactions = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/workspace/Maize project/MaizeReactionTable.txt");
for (my $i=1; $i < @{$maizereactions}; $i++) {
	my $items = [split(/\t/,$maizereactions->[$i])];
	$rxnec->{$items->[7]}->{$items->[1]} = 1;
}

my $reactionstbl = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/workspace/Maize project/Reactions.txt");
my $rxnhash = {};
my $keggpathways = {};
for (my $i=1; $i < @{$reactionstbl}; $i++) {
	my $items = [split(/\t/,$reactionstbl->[$i])];
	my $rxns = [split(/[,\s]+/,$items->[3])];
	for (my $j=0; $j < @{$rxns}; $j++) {
		if (!defined($rxnhash->{$rxns->[$j]})) {
			$rxnhash->{$rxns->[$j]} = {
				genes => {$items->[4] => 1},
				compartments => $items->[11],
				class => {$items->[0] => 1},
				subsystems => {$items->[1] => 1},
				roles => {$items->[2] => 1}
			};
		}
		$rxnhash->{$rxns->[$j]}->{class}->{$items->[0]} = 1;
		$rxnhash->{$rxns->[$j]}->{subsystems}->{$items->[1]} = 1;
		$rxnhash->{$rxns->[$j]}->{roles}->{$items->[2]} = 1;
		$rxnhash->{$rxns->[$j]}->{genes}->{$items->[4]} = 1;
		if (defined($rxndatahash->{$rxns->[$j]})) {
			if (defined($rxndatahash->{$rxns->[$j]}->{kegg_pathways})) {
				foreach my $mapid (@{$rxndatahash->{$rxns->[$j]}->{kegg_pathways}}) {
					$mapid =~ s/rn/map/;
					if (defined($kegghash->{$mapid})) {
						$keggpathways->{$kegghash->{$mapid}}->{$rxns->[$j]} = 1;
					}
				} 
			}
		}
	}
}

print "KEGG pathways\tKEGG ID\tEC number\tPlantSEED ID\tGenes\tCompartments\tEquation\tClass\tSubsystem\tRole\tComplete class\tHeterotoph class\tAutotroph class\n";
foreach my $pathway (keys(%{$keggpathways})) {
	foreach my $rxnid (keys(%{$keggpathways->{$pathway}})) {
		my $ec = "none";
		if (defined($rxndatahash->{$rxnid}->{ec_numbers})) {
			$ec = join("|",@{$rxndatahash->{$rxnid}->{ec_numbers}});
		} elsif (defined($rxnec->{$rxnid})) {
			$ec = join("|",keys(%{$rxnec->{$rxnid}}));
		}
		my $kegg = "none";
		if (defined($rxndatahash->{$rxnid}->{kegg_aliases})) {
			$kegg = join("|",@{$rxndatahash->{$rxnid}->{kegg_aliases}});
		}
		print $pathway."\t".
			$kegg."\t".
			$ec."\t".
			$rxnid."\t".
			join("|",keys(%{$rxnhash->{$rxnid}->{genes}}))."\t".
			$rxnhash->{$rxnid}->{compartments}."\t".
			$rxndatahash->{$rxnid}->{definition}."\t".
			join("|",keys(%{$rxnhash->{$rxnid}->{class}}))."\t".
			join("|",keys(%{$rxnhash->{$rxnid}->{subsystems}}))."\t".
			join("|",keys(%{$rxnhash->{$rxnid}->{roles}}))."\n";
	}
}