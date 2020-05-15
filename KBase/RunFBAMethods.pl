#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

#my $new_media_list = [
#	"filipeliu:narrative_1564175222344/Carbon-D-Glucose-iML1515",
#	"filipeliu:narrative_1564175222344/eae8f6da-cfd3-4b0b-bdcd-5f624918de7a-Carbon-D-Glucose-iML1515_v2",
#	"filipeliu:narrative_1564175222344/eae8f6da-cfd3-4b0b-bdcd-5f624918de7a-Carbon-D-Glucose-iML1515_v3"
#];
#
#my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
#	max_c_uptake => 60,
#	workspace => "filipeliu:narrative_1564175222344",
#	fbamodel_id => "iML1515.kb",
# 	fba_output_id => "iML1515.kb.mifba",
#	fbamodel_workspace => "filipeliu:narrative_1564175222344",
#	media_id_list => $new_media_list,
#	target_reaction => "bio1",
#	minimize_flux => 1
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

my $output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
	workspace => 49054,
	genome_id => "Lactococcus_lactis_Il1403.RAST", 
	fbamodel_output_id => "TestModel",
	media_id => "Carbon-D-Glucose",
	template_id => "auto",
	media_workspace => "KBaseMedia",
});

#my $output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 32329,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "New_keggro",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "classic",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["keggro"]
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 32329,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "New_SEED_keggro",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "classic",
#	anaerobe => 0,
#	use_annotated_functions => 1,
#	merge_all_annotations => 0,
#	source_ontology_list => ["keggro"]
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 32329,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "SSO_keggro",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["SSO","keggro"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_keggro",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["keggro"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_keggko",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["keggko"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_SSO_keggro",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["SSO","keggro"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_keggro_keggko",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["keggro","keggko"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_SSO_keggko",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["SSO","keggko"]
#});
#
#$output = Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
#	workspace => 30128,
#	genome_id => "PT32.4", 
#	fbamodel_output_id => "Merged_SSO_keggro_keggko",
#	media_id => "Carbon-D-Glucose",
#	template_id => "auto",
#	template_workspace => undef,
#	media_workspace => "KBaseMedia",
#	gapfill_model => 0,
#	mode => "new",
#	anaerobe => 0,
#	use_annotated_functions => 0,
#	merge_all_annotations => 0,
#	source_ontology_list => ["SSO","keggro","keggko"]
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

#my $output = $impl->export_fba_as_tsv_file({
#	input_ref => "43668/45"
#});

#my $modelid = "LuteusModel";
#my $expid = "VanillicAcid_SampleSet_TPM_ExpressionMatrix";
#my $mediaid = "Complete";
#my $mediaws = "KBaseMedia";
#my $ws = 39170;
#my $condition = "Mluteus_vanillicacid_1.fastq_reads_expression";

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

#my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
#	fbamodel_id => "Modified_Zmays_Model_Endosperm",
#	media_id => "KBaseMedia/Complete",
#	target_reaction => "bio1",
#	fba_output_id => "test",
#	fva => 0,
#	minimize_flux => 1,
#	thermodynamic_constraints => 0,
#	max_c_uptake => 60,
#	max_n_uptake => 90,
#	max_p_uptake => 15,
#	max_o_uptake => 40,
#	workspace => "Maize_Vitamin_Modeling"#,
#	#compute_characteristic_flux => 0,
#	#expseries_id => "Maize_qTeller",
#	#expseries_workspace => "Maize_Vitamin_Modeling",
#	#expression_condition => "Embryo_14DAP",
#	#characteristic_flux_file => "/Users/chenry/temp/fbajobs/eWZGbbEleC/MFAOutput/CharacteristicFluxes.txt"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_catalogue_all_loops_in_model({	
#	fbamodel_id => "PseudomonasSBW25.model",
#	workspace => "chenry:narrative_1515534936072",
#	fbamodel_output_id => "test"
#});

#my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
#	fbamodel_id => "CPCRW_v2",
#	workspace => "39182",
#	media_id => "RefGlucoseMinimal",
#	media_workspace => "39182",
#	target_reaction => "bio1",
#	fbamodel_output_id => "CPCRW_v2.1"
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
