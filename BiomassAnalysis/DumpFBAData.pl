use strict;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();

my $directory = "/Users/janakaanl/MSRepo/ResearchScripts/IGap/";
my $workspace = "janakakbase:narrative_1499012973666";
my $initial_models = Bio::KBase::kbaseenv::ws_client()->list_objects({
	workspaces => [$workspace],
	type => "KBaseFBA.FBAModel",
});
my $models;
for (my $i=0; $i < @{$initial_models}; $i++) {
	push(@{$models},$initial_models->[$i]->[1]);
}

my $initial_genomes = Bio::KBase::kbaseenv::ws_client()->list_objects({
       	workspaces => ["janakakbase:narrative_1499012973666"],
        type => "KBaseGenomes.Genome",
});

my $genomelist;
my $genome_hash = {};
for (my $i=0; $i < @{$initial_genomes}; $i++) {
        if ($initial_genomes->[$i]->[1] =~ m/(.+)\.RAST$/) {
                $genome_hash->{$initial_genomes->[$i]->[1]} = $1;
                push(@{$genomelist},$initial_genomes->[$i]->[1]);
        }
}

my $rxnhash;
my $mdlhash;
my $rxnobjhash;
for (my $i=0; $i < @{$models}; $i++) {
	print "Analyzing ".$models->[$i]."\n";
	my $model = $models->[$i];
	my $modelobj = $impl->util_store()->get_object($workspace."/".$model);
	my $rxncount = @{$modelobj->modelreactions()};
	my $cpdcount = @{$modelobj->modelcompounds()};
	my $genecount = @{$modelobj->features()};
	$mdlhash->{$models->[$i]} = {
		rxn => $rxncount,
		cpd => $cpdcount,
		genes => $genecount,
		gfrxn => 0,
		ssobj => 0,
		ssmedia => 0,
		mmactive => 0,
		comactive => 0,
		mmessential => 0,
		comessential => 0,
		mmflexible => 0,
		comflexible => 0,
		mmblocked => 0,
		comblocked => 0
	};
	my $rxns = $modelobj->modelreactions();
	for (my $j=0; $j < @{$rxns}; $j++) {
		$rxnobjhash->{$rxns->[$j]->id()} = $rxns->[$j];
		if ($rxns->[$j]->gprString() eq "Unknown") {
			$mdlhash->{$models->[$i]}->{gfrxn}++;
		}
		if (length($rxns->[$j]->gapfillString()) == 0) {
			$mdlhash->{$models->[$i]}->{gfrxn}++;
		}
		my $id = $rxns->[$j]->id();
		$id =~ s/_c0//;
		if (!defined($rxnhash->{$id})) {
			$rxnhash->{$id} = {
				equation => $rxns->[$j]->equation(),
				definition => $rxns->[$j]->definition(),
				name => $rxns->[$j]->name(),
				id => $id
			};
		}
		$rxnhash->{$id}->{$models->[$i]} = {
			gpr => $rxns->[$j]->gprString(),
			genelist => $rxns->[$j]->featureIDs(),
			direction => $rxns->[$j]->direction(),
			gapfillString => $rxns->[$j]->gapfillString()
		};
	}
	my $fbaworkspace = "janakakbase:narrative_1499012973666";
	my $fba = $models->[$i].".sensfba";
	my $fbaobj = $impl->util_store()->get_object($fbaworkspace."/".$fba);
	$mdlhash->{$models->[$i]}->{ssobj} = $fbaobj->objectiveValue();
	$mdlhash->{$models->[$i]}->{ssmedia} = $fbaobj->media()->_wsname();
	$rxns = $fbaobj->FBAReactionVariables();
	for (my $j=0; $j < @{$rxns}; $j++) {
		my $id = $rxns->[$j]->modelreaction()->id();
		$id =~ s/_c0//;
		$rxnhash->{$id}->{$models->[$i]}->{ssflux} = $rxns->[$j]->value();
		$rxnhash->{$id}->{$models->[$i]}->{biomassdep} = $rxns->[$j]->biomass_dependencies();
		$rxnhash->{$id}->{$models->[$i]}->{coupled_reactions} = $rxns->[$j]->coupled_reactions();
	}
	$fba = $models->[$i].".mmfva";
	$fbaobj = $impl->util_store()->get_object($fbaworkspace."/".$fba);
	$mdlhash->{$models->[$i]}->{mmobj} = $fbaobj->objectiveValue();
	$mdlhash->{$models->[$i]}->{mmmedia} = $fbaobj->media()->_wsname();
	my $cpds = $fbaobj->FBACompoundVariables();
	for (my $j=0; $j < @{$cpds}; $j++) {
		if ($cpds->[$j]->value() > 0.0000001 && $cpds->[$j]->modelcompound()->name() =~ m/glucose/i) {
			$mdlhash->{$models->[$i]}->{uptake}->{$cpds->[$j]->modelcompound()->name()} = $cpds->[$j]->value();
		} elsif ($cpds->[$j]->value() < -0.0000001 && $cpds->[$j]->modelcompound()->name() =~ m/lactate|formate|butyrate|acetate|succinate|propionate/i) {
			$mdlhash->{$models->[$i]}->{excretion}->{$cpds->[$j]->modelcompound()->name()} = $cpds->[$j]->value();
		}
	}
	$rxns = $fbaobj->FBAReactionVariables();
	for (my $j=0; $j < @{$rxns}; $j++) {
		my $id = $rxns->[$j]->modelreaction()->id();
		$id =~ s/_c0//;
		$rxnhash->{$id}->{$models->[$i]}->{mmflux} = $rxns->[$j]->value();
		$rxnhash->{$id}->{$models->[$i]}->{mmmax} = $rxns->[$j]->max();
		$rxnhash->{$id}->{$models->[$i]}->{mmmin} = $rxns->[$j]->min();
		if (abs($rxns->[$j]->value()) > 0.0000001) {
			$mdlhash->{$models->[$i]}->{mmacitve}++;
		}
		if ($rxns->[$j]->min() > 0.0000001) {
			$mdlhash->{$models->[$i]}->{mmessential}++;
			$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "p";
		} elsif ($rxns->[$j]->max() < -0.0000001) {
			$mdlhash->{$models->[$i]}->{mmessential}++;
			$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "n";
		} elsif ($rxns->[$j]->max() > 0.0000001) {
			if ($rxns->[$j]->min() < -0.0000001) {
				$mdlhash->{$models->[$i]}->{mmflexible}++;
				$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "v";
			} else {
				$mdlhash->{$models->[$i]}->{mmflexible}++;
				$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "pv";
			}
		} elsif ($rxns->[$j]->min() < -0.0000001) {
			$mdlhash->{$models->[$i]}->{mmflexible}++;
			$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "nv";
		} else {
			$mdlhash->{$models->[$i]}->{mmblocked}++;
			$rxnhash->{$id}->{$models->[$i]}->{mmclass} = "b";
		}
	}
	$fba = $models->[$i].".comfva";
	$fbaobj = $impl->util_store()->get_object($fbaworkspace."/".$fba);
	$mdlhash->{$models->[$i]}->{comobj} = $fbaobj->objectiveValue();
	$mdlhash->{$models->[$i]}->{commedia} = $fbaobj->media()->_wsname();
	$rxns = $fbaobj->FBAReactionVariables();
	for (my $j=0; $j < @{$rxns}; $j++) {
		my $id = $rxns->[$j]->modelreaction()->id();
		$id =~ s/_c0//;
		$rxnhash->{$id}->{$models->[$i]}->{comflux} = $rxns->[$j]->value();
		$rxnhash->{$id}->{$models->[$i]}->{commax} = $rxns->[$j]->max();
		$rxnhash->{$id}->{$models->[$i]}->{commin} = $rxns->[$j]->min();
		if (abs($rxns->[$j]->value()) > 0.0000001) {
			$mdlhash->{$genomelist->[$i]}->{comacitve}++;
			#$mdlhash->{$models->[$i]}->{comacitve}++;
		}
		if ($rxns->[$j]->min() > 0.0000001) {
			$mdlhash->{$models->[$i]}->{comessential}++;
			$rxnhash->{$id}->{$models->[$i]}->{comclass} = "p";
		} elsif ($rxns->[$j]->max() < -0.0000001) {
			$mdlhash->{$models->[$i]}->{comessential}++;
			$rxnhash->{$id}->{$models->[$i]}->{comclass} = "n";
		} elsif ($rxns->[$j]->max() > 0.0000001) {
			if ($rxns->[$j]->min() < -0.0000001) {
				$mdlhash->{$models->[$i]}->{comflexible}++;
				$rxnhash->{$id}->{$models->[$i]}->{comclass} = "v";
			} else {
				$mdlhash->{$models->[$i]}->{comflexible}++;
				$rxnhash->{$id}->{$models->[$i]}->{comclass} = "pv";
			}
		} elsif ($rxns->[$j]->min() < -0.0000001) {
			$mdlhash->{$models->[$i]}->{comflexible}++;
			$rxnhash->{$id}->{$models->[$i]}->{comclass} = "nv";
		} else {
			$mdlhash->{$models->[$i]}->{comblocked}++;
			$rxnhash->{$id}->{$models->[$i]}->{comclass} = "b";
		}
	}
}

my $json = Bio::KBase::ObjectAPI::utilities::TOJSON($rxnhash,1);
open(my $fout,">",$directory."ModelFBAData.json");
print $fout $json;
close($fout);
$json = Bio::KBase::ObjectAPI::utilities::TOJSON($mdlhash,1);
open(my $fout,">",$directory."MdlData.json");
print $fout $json;
close($fout);
