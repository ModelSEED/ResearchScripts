#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
use Bio::KBase::utilities;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $workspace = 51582;
my $growth = 0.5;
my $mediaid = "RichDefinedAnaerobic.media";
my $modelid = "Lactococcus_lactis_model";

my $media = $handler->util_get_object(Bio::KBase::utilities::buildref($mediaid,$workspace));
my $media_hash = {};
my $mediacpds = $media->mediacompounds();
foreach my $mediacpd (@{$mediacpds}) {
	$media_hash->{$mediacpd->compound()->id()} = $mediacpd;
}
delete $media_hash->{cpd00001};
delete $media_hash->{cpd00009};
delete $media_hash->{cpd00067};
delete $media_hash->{cpd00048};
delete $media_hash->{cpd00030};
delete $media_hash->{cpd00034};
delete $media_hash->{cpd00058};
delete $media_hash->{cpd00063};
delete $media_hash->{cpd00099};
delete $media_hash->{cpd00149};
delete $media_hash->{cpd00205};
delete $media_hash->{cpd00244};
delete $media_hash->{cpd00254};
delete $media_hash->{cpd00971};
delete $media_hash->{cpd10515};
delete $media_hash->{cpd10516};
delete $media_hash->{cpd10516};
delete $media_hash->{cpd00132};
delete $media_hash->{cpd00053};
delete $media_hash->{cpd00393};
delete $media_hash->{cpd00220};
delete $media_hash->{cpd00084};
delete $media_hash->{cpd00107};
delete $media_hash->{cpd00060};
delete $media_hash->{cpd00051};

my $model = $handler->util_get_object(Bio::KBase::utilities::buildref($modelid,$workspace));
my $biomasses = $model->biomasses();
foreach my $biomass (@{$biomasses}) {
	my $biocpds = $biomass->biomasscompounds();
	foreach my $biocpd (@{$biocpds}) {
		if (defined($media_hash->{$biocpd->modelcompound()->msid()})) {
			$media_hash->{$biocpd->modelcompound()->msid()}->maxFlux(-1*$biocpd->coefficient()*$growth);
		}
	}
}
$media_hash->{cpd00027}->maxFlux(20);

my $wsmeta = $handler->util_save_object($media,$workspace."/".$mediaid."2");

my $output = Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({	
	fbamodel_id => $modelid,
	media_id => $mediaid."2",
	target_reaction => "bio1",
	fba_output_id => "Lactis.fba",
	fva => 1,
	minimize_flux => 1,
	#max_c_uptake => 60,
	workspace => $workspace,
});

my $fluxdata = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/workspace/PlantSEED/CCRStudy/PreviousFluxes.txt");
my $fluxhash = {};
for (my $i=0; $i < @{$fluxdata}; $i++) {
	my $array = [split(/\t/,$fluxdata->[$i])];
	if (abs($array->[0]) > 0.0000001) {
		$fluxhash->{$array->[0]} = $array->[5];
	}
}

my $fba = $handler->util_get_object(Bio::KBase::utilities::buildref("Lactis.fba",$workspace));
my $rxns = $fba->FBAReactionVariables();
my $new = 0;
my $lost = 0;
my $newhash = {};
for (my $i=0; $i < @{$rxns}; $i++) {
	if (abs($rxns->[$i]->value()) > 0.0000001) {
		$newhash->{$rxns->[$i]->modelreaction()->id()} = $rxns->[$i]->value();
		if (!defined($fluxhash->{$rxns->[$i]->modelreaction()->id()})) {
			
		}
	}
}
foreach my $rxn (keys(%{$fluxhash})) {
	
}
