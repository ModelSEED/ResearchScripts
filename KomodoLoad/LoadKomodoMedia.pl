#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(parse_input_table save_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

my $cdata = parse_input_table($directory."/richcompounds.txt",[
	["STRMEDIAID",1],
	["STRMETABOLITEID",1],
	["STRSEMIUNIQUENAME",1],
	["STRCONCENTRATION",1],
	["STRUNITS",0],
	["DBLAMOUNT",0],
	["STRISRICH",0],
	["DBLTOTALMOLARWEIGHT",0],
	["STRSEEDMETABOLITES",0],
	["STRUNIQUENAME",0],
]);

my $metahash = {};
my $badmetabolites = {};
for (my $i=0; $i < @{$cdata}; $i++) {
	if (!defined($metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]})) {
		$metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]} = {
			id => $cdata->[$i]->[0],
			name => $cdata->[$i]->[2],
			concentration => $cdata->[$i]->[3],
			concentration_units => $cdata->[$i]->[4],
			molecular_weight => $cdata->[$i]->[7],
			is_rich => $cdata->[$i]->[6],
			associated_compounds => {}
		};
		if (!defined($cdata->[$i]->[8]) || length($cdata->[$i]->[8]) == 0) {
			$badmetabolites->{$cdata->[$i]->[9]}++;
		}
	}
}
foreach my $cpd (keys(%{$badmetabolites})) {
	print $cpd."\t".$badmetabolites->{$cpd}."\n";
}

$cdata = parse_input_table($directory."/mediacompounds.txt",[
	["STRMEDIAID",1],
	["STRMETABOLITEID",1],
	["STRSEEDMETABOLITE",1],
	["DBLMOLCOUNT",1],
	["DBLMOLPERLITER",0,0.001],
]);

my $mediacpdhash = {};
for (my $i=0; $i < @{$cdata}; $i++) {
	if (!defined($mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]})) {
		$mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]} = 0;
	}
	$mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]} += $cdata->[$i]->[4];
	$metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]}->{associated_compounds}->{$cdata->[$i]->[2]} = $cdata->[$i]->[3];
}

my $data = parse_input_table($directory."/media.txt",[
	["Media_Id",1],
	["Media_Name",1],
	["PH_info",0,"unknown"],
	["Atmosphere",0,"unknown"],
	["Atmosphere_additions",0,"unknown"],
	["rich_media",1],
	["taerobic",1],
	["protocol_link",1],
]);

my $idhash;
for (my $i=0; $i < @{$data}; $i++) {
	if ($data->[$i]->[5] eq "false") {
		$data->[$i]->[5] = 1;
	} else {
		$data->[$i]->[5] = 0;
	}
	if ($data->[$i]->[6] eq "false") {
		$data->[$i]->[6] = 0;
	} else {
		$data->[$i]->[6] = 1;
	}
	my $id = "k.".$data->[$i]->[0];
	my $array = [split(/_/,$id)];
	$id = $array->[0];
	if (defined($idhash->{$id})) {
		my $index = 1;
		while (defined($idhash->{$id.".".$index})) {
			$index++;
		}
		$id = $id.".".$index;
	}
	my $media = {
		id => $id,
		name => $data->[$i]->[1],
		source_id => $data->[$i]->[0],
		source => "DSMZ",
		protocol_link => $data->[$i]->[7],
		isDefined => $data->[$i]->[5],
		isMinimal => 0,
		isAerobic => $data->[$i]->[6],
		type => "Komodo",
		pH_data => $data->[$i]->[2],
		atmosphere => $data->[$i]->[3],
		atmosphere_addition => $data->[$i]->[4],
		reagents => [],
		mediacompounds => []
	};
	foreach my $reagent (keys(%{$metahash->{$data->[$i]->[0]}})) {
		push(@{$media->{reagents}},$metahash->{$data->[$i]->[0]}->{$reagent});
	}
	foreach my $reagent (keys(%{$mediacpdhash->{$data->[$i]->[0]}})) {
		push(@{$media->{mediacompounds}},{
			compound_ref => "921/29/8/compounds/id/".$reagent,
			concentration => $mediacpdhash->{$data->[$i]->[0]}->{$reagent},
			maxFlux => 100,
			minFlux => -100
		});
	}
	#save_workspace_object("KomodoMedia/".$id,$media,"KBaseBiochem.Media");
}		