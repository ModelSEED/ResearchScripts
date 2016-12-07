#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $fba = get_fba_client();
$fba->modify_reactions({
	model => "83333.1-testmodel",
	workspace => "chenrydemo",
	reactions => [[
		"rxn05292_c0",
		undef,
		"fig|83333.1.peg.1003",
		undef,
		undef,
		undef,
		undef
	]],
	output_id => "83333.1-testmodel.mod"
});