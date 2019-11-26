#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $genome_id = $ARGV[0];
my $workspace = $ARGV[1];

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $rxnhash = Bio::KBase::utilities::reaction_hash();

if ($genome_id =~ m/(.+)\.RAST$/) {
	my $name = $1;
	print "Now processing:".$genome_id."\n";
	my $datachannel = {};
	#Building base model which populates base ATP characteristics
	Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
		workspace => $workspace,
		genome_id => $genome_id,
		fbamodel_output_id => $name.".base",
		genome_workspace => $workspace,
		gapfill_model => 0,
	},$datachannel);
	my $auxo_output = Bio::KBase::ObjectAPI::functions::func_predict_auxotrophy_from_model({
		workspace => $workspace,
		fbamodel_id => $name.".base",
	},$datachannel);
	my $baseline_gapfilling = $datachannel->{fbamodel}->attributes()->{baseline_gapfilling};
	delete $datachannel->{fbamodel};
	my $gapfill_output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({
		workspace => $workspace,
		fbamodel_id => $name.".base",
		media_id => $name.".base".".auxo_media",
		fbamodel_output_id => $name.".gapfilled"
	},$datachannel);
	$datachannel->{fbamodel}->attributes()->{baseline_gapfilling} = $baseline_gapfilling;
	$datachannel->{fbamodel}->attributes()->{auxotrophy_gapfilling} = $gapfill_output->{number_gapfilled_reactions};
	$datachannel->{fbamodel}->attributes()->{gapfilled_atpprod} = $gapfill_output->{atpproduction};
	$datachannel->{fbamodel}->attributes()->{gapfilled_biomass} = $gapfill_output->{growth};
	my $fba_output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
		workspace => $workspace,
		fbamodel_id => $name.".gapfilled",
		fba_output_id => $name.".fba",
		media_id => $name.".base".".auxo_media",
		fva => 1,
		minimize_flux => 1,
		max_c_uptake => 30
	},$datachannel);	
	$datachannel->{fbamodel}->attributes()->{auxo_biomass} = $fba_output->{objective};
	my $fbaobj = $datachannel->{fba};
	my $rxnvar = $fbaobj->FBAReactionVariables();
	my $classhash = {};
	for (my $i=0; $i < @{$rxnvar}; $i++) {
		if ($rxnvar->[$i]->modelreaction_ref() =~ m/(rxn\d+)/) {
			$classhash->{$1}->{auxo} = $rxnvar->[$i]->{class};
		}
		if (!defined($datachannel->{fbamodel}->attributes()->{"auxo_class_".$rxnvar->[$i]->{class}})) {
			$datachannel->{fbamodel}->attributes()->{"auxo_class_".$rxnvar->[$i]->{class}} = 0;
		}
		$datachannel->{fbamodel}->attributes()->{"auxo_class_".$rxnvar->[$i]->{class}}++;
	}
	$fba_output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
		workspace => $workspace,
		fbamodel_id => $name.".gapfilled",
		fba_output_id => $name.".fba",
		fva => 1,
		minimize_flux => 1,
		max_c_uptake => 30
	},$datachannel);	
	$datachannel->{fbamodel}->attributes()->{complete_biomass} = $fba_output->{objective};
	my $fbaobj = $datachannel->{fba};
	my $rxnvar = $fbaobj->FBAReactionVariables();
	for (my $i=0; $i < @{$rxnvar}; $i++) {
		if ($rxnvar->[$i]->modelreaction_ref() =~ m/(rxn\d+)/) {
			$classhash->{$1}->{comp} = $rxnvar->[$i]->{class};
		}
		if (!defined($datachannel->{fbamodel}->attributes()->{"complete_class_".$rxnvar->[$i]->{class}})) {
			$datachannel->{fbamodel}->attributes()->{"complete_class_".$rxnvar->[$i]->{class}} = 0;
		}
		$datachannel->{fbamodel}->attributes()->{"complete_class_".$rxnvar->[$i]->{class}}++;
	}
	my $rxns = $datachannel->{fbamodel}->modelreactions();
	for (my $i=0; $i < @{$rxns}; $i++) {
		if ($rxns->[$i]->id() =~ m/(rxn\d+)/) {
			my $rxnid = $1;
			if (defined($rxnhash->{$rxnid}) && defined($rxnhash->{$rxnid}->{kegg_pathways})) {
				for (my $j=0; $j < @{$rxnhash->{$rxnid}->{kegg_pathways}}; $j++) {
					if (!defined($datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]})) {
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]} = {
							rxn => 0, nonblocked => 0, gf => 0,rxns => {}
						};
					}
					if (length($rxns->[$i]->gapfillString()) > 0) {
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{gf}++;
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{rxns}->{$rxnid} = "g";
					} else {
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{rxn}++;
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{rxns}->{$rxnid} = "n";
					}
					if (defined($classhash->{$rxnid}->{comp}) && $classhash->{$rxnid}->{comp} ne "Blocked") {
						$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{nonblocked}++;
						if ($datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{rxns}->{$rxnid} ne "g") {
							$datachannel->{fbamodel}->attributes()->{"pathways_".$rxnhash->{$rxnid}->{kegg_pathways}->[$j]}->{rxns}->{$rxnid} = "a";
						}
					}
				}
			}
		}
	}
	foreach my $cpd (keys(%{$auxo_output->{auxotrophy_data}})) {
		$datachannel->{fbamodel}->attributes()->{"auxotrophy_".$cpd} = $auxo_output->{auxotrophy_data}->{$cpd}->{name}.":".
			$auxo_output->{auxotrophy_data}->{$cpd}->{totalrxn}.":".
			$auxo_output->{auxotrophy_data}->{$cpd}->{gfrxn}.":".
			$auxo_output->{auxotrophy_data}->{$cpd}->{auxotrophic};
	}
	my $wsmeta = $impl->util_save_object($datachannel->{fbamodel},Bio::KBase::utilities::buildref($name.".gapfilled",$workspace));
}