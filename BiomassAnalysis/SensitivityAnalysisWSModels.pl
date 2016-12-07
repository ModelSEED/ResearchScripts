#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];
my $media = $ARGV[1];
my $procindex = $ARGV[2];
my $numprocs = $ARGV[3];
my $fba = get_fba_client();
my $models = get_ws_objects_list($workspace,"KBaseFBA.FBAModel");
for (my $i=0; $i < @{$models}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		print $models->[$i]->[1]."\n";
		my $output = $fba->reaction_sensitivity_analysis({
			model => $models->[$i]->[1],
			workspace => $workspace,
			integrate_solution => 1,
			rxnsens_uid => $models->[$i]->[1].".rs",
			media => $media,
			media_ws => "KBaseMedia",
			delete_essential_reactions => 1,
			objective_sensitivity_only => 1,
			solver => "cplex"
		});
		printObjectInfo($output);
	}
}