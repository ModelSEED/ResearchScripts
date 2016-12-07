#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

#my $prefix = "OSCAR";

my $ws = "chenry:1431835409789";
(my $data,my $prov) = get_workspace_object($ws."/iBsuDSM_merged");
#Adding cpd16336
push(@{$data->{modelcompounds}},{
	id => "cpd16336_c0",
	name => "isoprene_c0",
	compound_ref => "489/6/1/compounds/id/cpd16336",
	modelcompartment_ref => "~/modelcompartments/id/c0",
	charge => 0,
	formula => "C5H8"
});
push(@{$data->{modelcompounds}},{
	id => "cpd16336_e0",
	name => "isoprene_e0",
	compound_ref => "489/6/1/compounds/id/cpd16336",
	modelcompartment_ref => "~/modelcompartments/id/e0",
	charge => 0,
	formula => "C5H8"
});
#Adding rxn11733
push(@{$data->{modelreactions}},{
	id => "rxn11733_c0",
	name => "dimethylallyl-diphosphate diphosphate-lyase",
	pathway => "Isoprene biosynthesis and excretion",
	reaction_ref => "489/6/1/reactions/id/rxn11733",
	direction => ">",
	protons => 0,
	probability => 1,
	modelcompartment_ref => "~/modelcompartments/id/c0",
	modelReactionReagents => [
		{"coefficient" => -1,"modelcompound_ref" => "~/modelcompounds/id/cpd00202_c0"},
		{"coefficient" => 1,"modelcompound_ref" => "~/modelcompounds/id/cpd00012_c0"},
		{"coefficient" => 1,"modelcompound_ref" => "~/modelcompounds/id/cpd16336_c0"}
	],
	modelReactionProteins => []
});
#Adding rxn13765
push(@{$data->{modelreactions}},{
	id => "rxn13765_c0",
	name => "Isoprene diffusion",
	pathway => "Isoprene biosynthesis and excretion",
	reaction_ref => "489/6/1/reactions/id/rxn13765",
	direction => "=",
	protons => 0,
	probability => 1,
	modelcompartment_ref => "~/modelcompartments/id/c0",
	modelReactionReagents => [
		{"coefficient" => -1,"modelcompound_ref" => "~/modelcompounds/id/cpd16336_c0"},
		{"coefficient" => 1,"modelcompound_ref" => "~/modelcompounds/id/cpd16336_e0"}
	],
	modelReactionProteins => []
});
save_workspace_object($ws."/iBsuDSM_isoprene",$data,"KBaseFBA.FBAModel");