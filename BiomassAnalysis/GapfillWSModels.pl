#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];
my $procindex = $ARGV[1];
my $numprocs = $ARGV[2];
my $models = get_ws_objects_list($workspace,"KBaseFBA.FBAModel");
for (my $i=0; $i < @{$models}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		my $fba = get_fba_client();
		my $output = $fba->adjust_biomass_reaction({
			compounds => ["cpd11715","cpd11746","cpd09680"],
			model => $models->[$i]->[1],
			workspace => "chenry:BiomassAnalysisMMGF",
			coefficients => [0,0,0],
    		compartments => ["c","c","c"],
    		compartmentIndecies => [0,0,0],
    	});
    	printObjectInfo($output);
		$output = $fba->gapfill_model({
			model => $models->[$i]->[1],
			model_workspace => $workspace,
			workspace => "chenry:BiomassAnalysisMMGF",
			integrate_solution => 1,
			out_model => $models->[$i]->[1].".gf",
			solver => "CPLEX",
			fastgapfill => 1,
			formulation => {
				formulation => {
					media => "Carbon-D-Glucose",
					media_workspace => "KBaseMedia"
				}
			}
		});
		printObjectInfo($output);
		$output = $fba->gapfill_model({
			model => $models->[$i]->[1],
			model_workspace => $workspace,
			workspace => "chenry:BiomassAnalysisCMGF",
			integrate_solution => 1,
			out_model => $models->[$i]->[1].".gf",
			solver => "CPLEX",
			fastgapfill => 1,
		});
		printObjectInfo($output);
	}
}