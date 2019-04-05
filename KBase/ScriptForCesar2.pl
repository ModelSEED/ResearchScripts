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
my $ws = Bio::KBase::kbaseenv::ws_client();

##CELL 2019-04-02 CODE FOR CHRIS
my $my_ws = "cesarcardonau:narrative_1553799037978";

#RUNNING THE CODE THRES 09 #04-02-2019
##Pre:Requires objects "CM_thr_09" and "AnyC.PS" and "AnyN.PS"
#my $media_C = &run_gf_for_PS({
#	"modelid" => "CM_thr_09",
#	"PS_name" => "AnyC.PS"
#});
#my $media_N = &run_gf_for_PS({
#	"modelid" => "CM_thr_09",
#	"PS_name" => "AnyN.PS"
#});
#&run_gf({
#	modelid => "CM_thr_09",
#	media => "KBaseMedia/ArgonneLBMedia",
#	modelid_out_suffix => "_gf_ArgonneLBMedia"
#});

##ADDITIONAL CODE THRES 07 #04-03-2019
##Pre:Requires objects "SBacteria_thr_07" and "AnyC.PS" and "AnyN.PS"
&run_gf_complete({
	modelid => "SBacteria_thr_07"
});
&run_gf_g10({
	modelid => "SBacteria_thr_07_gf_complete"
});
$ws->copy_object({
	from => {"ref" => $my_ws.'/SBacteria_thr_07_gf_complete_gf_g10'},
	to => {'ref' => $my_ws."/CM_thr_07"}
});
my $media_C = &run_gf_for_PS({
	"modelid" => "CM_thr_07",
	"PS_name" => "AnyC.PS"
});
my $media_N = &run_gf_for_PS({
	"modelid" => "CM_thr_07",
	"PS_name" => "AnyN.PS"
});
&run_gf({
	modelid => "CM_thr_07",
	media => "KBaseMedia/ArgonneLBMedia",
	modelid_out_suffix => "_gf_ArgonneLBMedia"
});

##ADDITIONAL CODE THRES 05 #04-03-2019
##Pre:Requires objects "SBacteria_thr_05" and "AnyC.PS" and "AnyN.PS"
#&run_gf_complete({
#	modelid => "SBacteria_thr_05"
#});
#&run_gf_g10({
#	modelid => "SBacteria_thr_05_gf_complete"
#});
#$ws->copy_object({
#	from => {"ref" => $my_ws.'/SBacteria_thr_05_gf_complete_gf_g10'},
#	to => {'ref' => $my_ws."/CM_thr_05"}
#});
#my $media_C = &run_gf_for_PS({
#	"modelid" => "CM_thr_05",
#	"PS_name" => "AnyC.PS"
#});
#my $media_N = &run_gf_for_PS({
#	"modelid" => "CM_thr_05",
#	"PS_name" => "AnyN.PS"
#});
#&run_gf({
#	modelid => "CM_thr_05",
#	media => "KBaseMedia/ArgonneLBMedia",
#	modelid_out_suffix => "_gf_ArgonneLBMedia"
#});

sub run_gf {
	my ($params) = @_;
	$params = Bio::KBase::utilities::args($params,["modelid","media"],{
		modelid_out_suffix => '_gf_complete'
	});
	my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
		fbamodel_id => $params->{modelid},
		workspace => $my_ws,
		media_id => $params->{media},
		target_reaction => "bio1",
		fbamodel_output_id => $params->{modelid}.$params->{modelid_out_suffix}
	});
	return $params->{modelid}.$params->{modelid_out_suffix};
}

sub run_gf_g10 { 
    my ($params) = @_;
	$params = Bio::KBase::utilities::args($params,["modelid"],{
		modelid_out_suffix => '_gf_g10'
	});
	$params->{media} = "C-D-Glucose10";
    &run_gf($params);
}

sub run_gf_complete {
	my ($params) = @_;
	$params = Bio::KBase::utilities::args($params,["modelid"],{
		modelid_out_suffix => '_gf_complete'
	});
	$params->{media} = "KBaseMedia/Complete";
    &run_gf($params);
}

sub run_gf_for_PS {
	my ($params) = @_;
	$params = Bio::KBase::utilities::args($params,["modelid","PS_name"],{
		names_only => 0
	});
    my $models_run = Bio::KBase::kbaseenv::list_objects({'workspaces' => [$my_ws], 'type' => 'KBaseFBA.FBAModel'});
    my $PS_obj = Bio::KBase::kbaseenv::get_object($my_ws,$params->{PS_name});
	my $media_list_KB = Bio::KBase::kbaseenv::list_objects({'workspaces' => ['KBaseMedia'], 'type' => 'KBaseBiochem.Media'});
    my $media_found = [];
    #loop through all our media in the phenotype set (PS)
    for (my $i=0; $i < @{$PS_obj->{'phenotypes'}}; $i++) {
        my $mediaPS_already_run = 0;
        my $mediaPS_i=$PS_obj->{"phenotypes"}->[$i];
        foreach my $mediaKbase (@{$media_list_KB}) {
            my $array = [split(/\//,$mediaPS_i->{"media_ref"})];
            if ($mediaKbase->[0] eq $array->[1]) {
                my $media_name = $mediaKbase->[1];
                push(@{$media_found},$media_name);
                my $gf_model_name = $params->{modelid}."_gf_".$media_name;
                foreach my $model_run_i (@{$models_run}) {
                    if ($model_run_i->[1] eq $gf_model_name) {
                        $mediaPS_already_run = 1;
                        last;
                    }
                }
                if ($mediaPS_already_run == 0 && $params->{names_only} == 0) {
                    &run_gf({
						modelid => $params->{modelid},
						media => 'KBaseMedia/'.$media_name,
						modelid_out_suffix => "_gf_".$media_name
					});
                }
            }
        }
    }
    return $media_found;
}