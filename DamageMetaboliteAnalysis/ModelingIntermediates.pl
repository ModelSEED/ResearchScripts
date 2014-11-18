#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $procindex = $ARGV[0];
my $numprocs = $ARGV[1];

my $fba = get_fba_client();
my $fbaout;
my $media = get_ws_objects_list("KBaseMedia","KBaseBiochem.Media");
for (my $i=0; $i < @{$media}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		my $output = $fba->runfba({
			model => "iBsu1103",
			model_workspace => "jplfaria:modelingtranscriptomics",
			workspace => "chenry:MetabliteModelingAnalysis",
			fva => 1,
			formulation => {
				media => $media->[$i]->[1],
				media_workspace => "KBaseMedia",
			}
	   	});
	   	if ($output->[10]->{Objective} > 0.0000001) {
	   		push(@{$fbaout->{iBsu1103}},$output);
	   	}  	
	}
}
for (my $i=0; $i < @{$media}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		my $output = $fba->runfba({
			model => "83333.1.fbamdl",
			model_workspace => "chenrydemo",
			workspace => "chenry:MetabliteModelingAnalysis",
			fva => 1,
			formulation => {
				media => $media->[$i]->[1],
				media_workspace => "KBaseMedia",
			}
	   	});
	   	if ($output->[10]->{Objective} > 0.0000001) {
	   		push(@{$fbaout->{"83333.1.fbamdl"}},$output);
	   	}  	
	}
}		