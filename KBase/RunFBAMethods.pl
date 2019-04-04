#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

#my $modelid = "LuteusModel";
#my $expid = "VanillicAcid_SampleSet_TPM_ExpressionMatrix";
#my $mediaid = "Complete";
#my $mediaws = "KBaseMedia";
#my $ws = 39170;
#my $condition = "Mluteus_vanillicacid_1.fastq_reads_expression";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

#my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
#	fba_output_id => $modelid.".fba",
#	fbamodel_id => $modelid,
#	workspace => $ws,
#	media_id => $mediaid,
#	media_workspace => "KBaseMedia",
#	target_reaction => "bio1",
#	expseries_id => $expid,
#	expression_condition => $condition,
#	exp_threshold_percentile => 0.5,
#	exp_threshold_margin => 0.1,
#	activation_coefficient => 0.5,
#	omega => 0,
#	objective_fraction => 0.1,
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_build_metagenome_metabolic_model({
#	workspace => 24430,
#	input_ref => "SoilSFA_WA_nrKO.genome.assembly",
#	fbamodel_output_id => "test_metagenome_model",
#	media_id => "Carbon-D-Glucose",
#	media_workspace => "KBaseMedia",
#	gapfill_model => 1,
#	gff_file => "/Users/chenry/workspace/PNNLSFA/SoilSFA_WA_nrKO.gff"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_fit_exometabolite_data({
#	workspace => 39597,
#	fbamodel_id => "P_flourescens_GW456-L13.mdl",
#	media_ref => "KBaseMedia/Carbon-D-Glucose",
#	metabolite_condition => "R2A_Pseudomonas_GW456-L13",
#	metabolite_matrix_ref => "39477/16/1",
#	target_reaction => "bio1",
#	fbamodel_output_id => "P_flourescens_GW456-L13.recon",
#	atp_production_check => 1
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
#	fbamodel_id => "cesarcardonau:narrative_1550613778929/Sloan.model_gf_Carbon-a-D-Glucose",
#	media_id => "KBaseMedia/Carbon-a-D-Glucose",
#	target_reaction => "bio1",
#	fba_output_id => "test",
#	fva => 0,
#	minimize_flux => 1,
#	predict_community_composition => 1,
#	max_c_uptake => 60,
#	max_n_uptake => 100,
#	workspace => "cesarcardonau:narrative_1550613778929"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
#	fbamodel_id => "PseudomonasSBW25.model",
#	media_id => "KBaseMedia/Complete",
#	target_reaction => "bio1",
#	fba_output_id => "test",
#	fva => 0,
#	minimize_flux => 1,
#	thermodynamic_constraints => 1,
#	max_c_uptake => 60,
#	workspace => "chenry:narrative_1515534936072"
#});

my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
	fbamodel_id => "Modified_Zmays_Model_Endosperm",
	media_id => "KBaseMedia/Complete",
	target_reaction => "bio1",
	fba_output_id => "test",
	fva => 0,
	minimize_flux => 1,
	thermodynamic_constraints => 0,
	max_c_uptake => 60,
	max_n_uptake => 90,
	max_p_uptake => 15,
	max_o_uptake => 40,
	workspace => "Maize_Vitamin_Modeling"#,
	#compute_characteristic_flux => 0,
	#expseries_id => "Maize_qTeller",
	#expseries_workspace => "Maize_Vitamin_Modeling",
	#expression_condition => "Embryo_14DAP",
	#characteristic_flux_file => "/Users/chenry/temp/fbajobs/eWZGbbEleC/MFAOutput/CharacteristicFluxes.txt"
});

#my $output = Bio::KBase::ObjectAPI::functions::func_catalogue_all_loops_in_model({	
#	fbamodel_id => "PseudomonasSBW25.model",
#	workspace => "chenry:narrative_1515534936072",
#	fbamodel_output_id => "test"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
#	fbamodel_id => "Sloan_050.model",
#	workspace => "40177",
#	media_id => "C-D-Glucose10",
#	media_workspace => "40177",
#	target_reaction => "bio1",
#	fbamodel_output_id => "Sloan_050.model.gf"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_predict_metabolite_biosynthesis_pathway({
#	target_metabolite_list => ["cpd00703_c0"],
#	fba_output_id => "SBW25_pathway_cpd00703",
#	fbamodel_id => "PseudomonasSBW25.model",
#	workspace => "chenry:narrative_1515534936072",
#	media_id => "Complete",
#	media_workspace => "KBaseMedia",
#	feature_ko_list => [],
#	reaction_ko_list => [],
#	custom_bound_list => [],
#	media_supplement_list => [],
#	source_metabolite_list => ["cpd00065_c0"],
#});
