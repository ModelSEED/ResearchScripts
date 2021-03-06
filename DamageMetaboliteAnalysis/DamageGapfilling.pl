#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw( fbaws get_fba_client runFBACommand );

my $fba = get_fba_client();

my $output = $fba->gapfill_model({
	model => "Reaction_Pair_RXN0001_c0_RXN0002_c0_GPP",
	workspace => "janakakbase:MetRepair_Pairs",
	integrate_solution => 1,
	model_workspace => "janakakbase:MetRepair_Pairs",
	out_model => "Reaction_Pair_RXN0001_c0_RXN0002_c0_GPP.gf",
	target_reactions => ["RXN0001_c0"],
	solver => "cplex",
	fastgapfill => 1,
	source_model => "AllDR_Jan9_GP",
	source_model_ws => "janakakbase:MetRepair"
});

my $output = $fba->runfba({
	model => "Reaction_Pair_RXN0001_c0_RXN0002_c0_GPP.gf",
	workspace => "janakakbase:MetRepair_Pairs",
	formulation => {
		bounds => [[0.2,0.2,"flux","RXN0001_c0"]]
	},
	solver => "cplex",
	fba => "Reaction_Pair_RXN0001_c0_RXN0002_c0_GPP.gf.fba",
	target_reactions => ["RXN0001_c0"],
	solver => "cplex",
	fastgapfill => 1,
	source_model => "AllDR_Jan9_GP",
	source_model_ws => "janakakbase:MetRepair"
});