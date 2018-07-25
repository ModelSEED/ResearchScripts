#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
use Bio::KBase::utilities;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $directory = $ARGV[0];
if (!defined($directory) || $directory eq "") {
	$directory = "/Users/chenry/Dropbox/workspace/Cesar/";
}
my $modelfile = $directory."patric-16S-selected-blasted-filtered-annotations.map";
my $workspace = "chenry:narrative_1530974574420";
my $mdldata = Bio::KBase::ObjectAPI::utilities::LOADFILE($modelfile);
my $mdlhash;
for (my $i=0; $i < @{$mdldata}; $i++) {
	my $array = [split(/,/,$mdldata->[$i])];
	$mdlhash->{$array->[0]} = 1;
}

my $params = {
	workspace => "chenry:narrative_1530974574420",
	fbamodel_id => undef,
	fba_output_id => undef,
	media_id => "Carbon-D-Glucose",
	media_workspace => "KBaseMedia"
};
foreach my $mdlid (keys(%{$mdlhash})) {
	if (!-e $directory."FBA/".$mdlid.".touch") {
		$params->{fbamodel_id} = $mdlid;
		$params->{fba_output_id} = $mdlid.".metabolite_interaction.fba";
		$params = Bio::KBase::utilities::args($params,["workspace","fbamodel_id"],{
			fbamodel_workspace => $params->{workspace},
			media_id => undef,
			media_workspace => $params->{workspace},
			probanno_id => undef,
			probanno_workspace => $params->{workspace},
			target_reaction => "bio1",
			fbamodel_output_id => $params->{fbamodel_id},
			thermodynamic_constraints => 0,
			comprehensive_gapfill => 0,
			source_fbamodel_id => undef,
			source_fbamodel_workspace => $params->{workspace},
			feature_ko_list => [],
			reaction_ko_list => [],
			custom_bound_list => [],
			media_supplement_list => [],
			expseries_id => undef,
			expseries_workspace => $params->{workspace},
			expression_condition => undef,
			exp_threshold_percentile => 0.5,
			exp_threshold_margin => 0.1,
			activation_coefficient => 0.5,
			omega => 0,
			objective_fraction => 0,
			minimum_target_flux => 0.1,
			number_of_solutions => 1,
			gapfill_output_id => undef
		});
		my $printreport = 1;
		my $htmlreport = "";
		if (defined($params->{reaction_ko_list}) && ref($params->{reaction_ko_list}) ne "ARRAY") {
			if (length($params->{reaction_ko_list}) > 0) {
				$params->{reaction_ko_list} = [split(/,/,$params->{reaction_ko_list})];
			} else {
				 $params->{reaction_ko_list} = [];
			}
		}
		my $model;
		if (!defined($model)) {
			$handler->util_log("Retrieving model.");
			$model = $handler->util_get_object(Bio::KBase::utilities::buildref($params->{fbamodel_id},$params->{fbamodel_workspace}));
			$htmlreport .= Bio::KBase::utilities::style()."<div style=\"height: 200px; overflow-y: scroll;\"><p>The genome-scale metabolic model ".$params->{fbamodel_id}." was gapfilled";
		} else {
			$printreport = 0;
			$htmlreport .= "<p>The model ".$params->{fbamodel_id}." was gapfilled";
		}
		if (!defined($params->{media_id})) {
			if ($model->genome()->domain() eq "Plant" || $model->genome()->taxonomy() =~ /viridiplantae/i) {
				$params->{media_id} = Bio::KBase::utilities::conf("ModelSEED","default_plant_media");
			} else {
				$params->{default_max_uptake} = 100;
				$params->{media_id} = Bio::KBase::utilities::conf("ModelSEED","default_microbial_media");
			}
			$params->{media_workspace} = Bio::KBase::utilities::conf("ModelSEED","default_media_workspace");
		}
		$htmlreport .= " in ".$params->{media_id}." media to force a minimum flux of ".$params->{minimum_target_flux}." through the ".$params->{target_reaction}." reaction.";
		$handler->util_log("Retrieving ".$params->{media_id}." media.");
		my $media = $handler->util_get_object(Bio::KBase::utilities::buildref($params->{media_id},$params->{media_workspace}));
		$handler->util_log("Preparing flux balance analysis problem.");
		my $source_model;
		if (defined($params->{source_fbamodel_id}) && !defined($source_model)) {
			$htmlreport .= " During the gapfilling, the source biochemistry database was augmented with all the reactions contained in the existing ".$params->{source_fbamodel_id}." model.";
			$source_model = $handler->util_get_object(Bio::KBase::utilities::buildref($params->{source_fbamodel_id},$params->{source_fbamodel_workspace}));	
		}
		my $gfs = $model->gapfillings();
		my $currentid = 0;
		for (my $j=0; $j < @{$gfs}; $j++) {
			if ($gfs->[$j]->id() =~ m/gf\.(\d+)$/) {
				if ($1 >= $currentid) {
					$currentid = $1+1;
				}
			}
		}
		my $gfid = "gf.".$currentid;
		my $fba = Bio::KBase::ObjectAPI::functions::util_build_fba($params,$model,$media,$params->{fbamodel_output_id}.".".$gfid,1,1,$source_model,1);
		$fba->parameters()->{"Metabolite interaction analysis"} = 1;
		$fba->parameters()->{"Target metabolite list"} = "cpd01682_c0;cpd00128_c0;cpd00540_c0;cpd01220_c0;cpd19247_c0;cpd00098_c0;cpd09874_c0;cpd26883_c0;cpd06619_c0;cpd23360_c0;cpd00089_c0;cpd00053_c0;cpd00119_c0;cpd00398_c0;cpd00133_c0;cpd00218_c0;cpd03161_c0;cpd00323_c0;cpd01209_c0;cpd10138_c0;cpd00737_c0;cpd00051_c0;cpd00586_c0;cpd07061_c0;cpd19072_c0;cpd19822_c0;cpd02601_c0";
		$fba->parameters()->{"Core metabolite list"} = "cpd00020_c0;cpd00024_c0;cpd00169_c0;cpd00102_c0;cpd00072_c0;cpd00032_c0;cpd00079_c0;cpd00022_c0;cpd00236_c0;cpd00101_c0;cpd00061_c0";
		$fba->fluxMinimization(1);
		$handler->util_log("Running flux balance analysis problem.");
		my $objective = $fba->runFBA();
		$fba->toJSON({pp => 1});
		$handler->util_log("Saving FBA results.");
		$fba->id($params->{fba_output_id});
		my $wsmeta = $handler->util_save_object($fba,Bio::KBase::utilities::buildref($params->{fba_output_id},$params->{workspace}),{type => "KBaseFBA.FBA"});
		system("touch ".$directory."FBA/".$mdlid.".touch");
		exit();
	}
}