#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $model = $ARGV[0];
my $fba = get_fba_client();
my $media = get_ws_objects_list("RhodoMedia","KBaseBiochem.Media");
for (my $i=0; $i < @{$media}; $i++) {
	my $output = $fba->runfba({
		model => $model,
		model_workspace => "MetaboliteEssentialityAnalysis",
		workspace => "MetaboliteEssentialityAnalysis",
		fva => 1,
		formulation => {
			media => $media->[$i]->[1],
			media_workspace => "RhodoMedia",
		}
   	});
   	if ($output->[10]->{Objective} > 0.0000001) {
   		print join("\t",@{$output})."\n";
   	}
}