#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(parse_input_table save_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

my $fba = get_fba_client();
my $output = $fba->get_biochemistry({});
my $cpds = $output->{compounds};
my $compounds = $fba->get_compounds({
	compounds => $cpds
});
my $cpdhash = {};
for (my $i=0; $i < @{$compounds}; $i++) {
	$cpdhash->{$compounds->[$i]->{id}} = $compounds->[$i]->{formula};
}

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

my $subs = [
"rich-seawater#	cpd00001;cpd00099;cpd00971",
"H3BO#	cpd00067|3",
"H3BO4#	cpd00067|3",
"C12H22O11#	cpd00076",
"SnCl2#	cpd00099",
"rich-DNA#	cpd00115;cpd00357;cpd00241;cpd00356",
"phenol-solution#	cpd00127",
"(NH4)6 Mo7O24#	cpd00131|7;cpd00013|6",
"rich-charcoal#	cpd00242",
"polygalacturonate#	cpd00280",
"rich-casamino#	cpd00322;cpd00156;cpd00069;cpd00065;cpd00161;cpd00054;cpd00129;cpd00066;cpd00017;cpd00039;cpd00107;cpd00119;cpd00023;cpd00053;cpd00084;cpd00041;cpd00132;cpd00035;cpd00051;cpd00033",
"rich-meat#	cpd00322;cpd00156;cpd00069;cpd00065;cpd00161;cpd00054;cpd00129;cpd00066;cpd00017;cpd00039;cpd00107;cpd00119;cpd00023;cpd00053;cpd00084;cpd00041;cpd00132;cpd00035;cpd00051;cpd00033",
"Streptomycin#	cpd00328",
"Naphthalene#	cpd00618",
"6-methylnicotinate#	cpd00737",
"3-hydroxybutyrate#	cpd00797",
"NaVO3#	cpd00971",
"2;4-dihydroxybenzoate-Na salt#	cpd00971",
"Na2MoO2#	cpd00971|2",
"o-nitrophenol#	cpd01361",
"haemin#	cpd01476",
"1;4-Naphthaquinone#	cpd01703",
"2#rich-peptone#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"rich-serum#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"rich-peptone#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"rich-blood#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"rich-plant#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"rich-rumen#	cpd03424;cpd00215;cpd00028;cpd10515;cpd00030;cpd00149;cpd00058;cpd00099;cpd00007;cpd00034;cpd00156;cpd00249;cpd00092;cpd00069;cpd00065;cpd00184;cpd00161;cpd00048;cpd00054;cpd00220;cpd00129;cpd00644;cpd00009;cpd00066;cpd00218;cpd00971;cpd00254;cpd00060;cpd00039;cpd00107;cpd00205;cpd00246;cpd00322;cpd00226;cpd00119;cpd00531;cpd00001;cpd00067;cpd00033;cpd00023;cpd00027;cpd00393;cpd10516;cpd00654;cpd00438;cpd00381;cpd11595;cpd01012;cpd00063;cpd00041;cpd01048;cpd00051;cpd00035;cpd00182;cpd00311;cpd00126;cpd00018;cpd00091;cpd00046;cpd00793;cpd00541;cpd00239",
"Antipyrine#	cpd04709",
"ACES#	cpd11966",
"alum#	cpd12855",
"aluminum chloride#	cpd12855;cpd00099"
];

my $universals = [qw(
cpd00011
cpd10516
cpd00067
cpd00001
cpd00205
cpd00030
cpd00034
cpd10515
)];

my $methash;
for (my $i=0; $i < @{$subs}; $i++) {
	my $array = [split(/\t/,$subs->[$i])];
	my $cpdarray = [split(/;/,$array->[1])];
	for (my $j=0; $j < @{$cpdarray}; $j++) {
		my $newarray = [split(/\|/,$cpdarray->[$j])];
		if (defined($newarray->[1])) {
			$methash->{$array->[0]}->{$newarray->[0]} = $newarray->[1];
		} else {
			$methash->{$array->[0]}->{$newarray->[0]} = 1;
		}	
	}
}

my $metahash = {};
my $badmetabolites = {};
my $mediacpdhash = {};
for (my $i=0; $i < @{$cdata}; $i++) {
	if (!defined($metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]})) {
		if ($cdata->[$i]->[6] eq "FALSE") {
			$cdata->[$i]->[6] = 0;
		} else {
			$cdata->[$i]->[6] = 1;
		}
		if ($cdata->[$i]->[3] =~ m/(\d+\.*\d*)/) {
			$cdata->[$i]->[3] = $1;
		} elsif ($cdata->[$i]->[3] =~ m/trace/) {
			$cdata->[$i]->[3] = 0.00001;
		} elsif ($cdata->[$i]->[3] =~ m/substr/) {
			$cdata->[$i]->[3] = 1;
		}
		$metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]} = {
			id => $cdata->[$i]->[0],
			name => $cdata->[$i]->[2],
			concentration => $cdata->[$i]->[3]+0,
			concentration_units => $cdata->[$i]->[4],
			molecular_weight => $cdata->[$i]->[7]+0,
			is_rich => $cdata->[$i]->[6]+0,
			associated_compounds => {}
		};
		if (!defined($cdata->[$i]->[8]) || length($cdata->[$i]->[8]) == 0) {
			if (defined($methash->{$cdata->[$i]->[9]})) {
				foreach my $cpd (keys(%{$methash->{$cdata->[$i]->[9]}})) {
					$metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]}->{associated_compounds}->{$cpd} = $methash->{$cdata->[$i]->[9]}->{$cpd}+0;
					if (!defined($mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]})) {
						$mediacpdhash->{$cdata->[$i]->[0]}->{$cpd} = 0;
					}
					$mediacpdhash->{$cdata->[$i]->[0]}->{$cpd} += $methash->{$cdata->[$i]->[9]}->{$cpd}*$cdata->[$i]->[3];
				}
			}
		}
	}
}

$cdata = parse_input_table($directory."/mediacompounds.txt",[
	["STRMEDIAID",1],
	["STRMETABOLITEID",1],
	["STRSEEDMETABOLITE",1],
	["DBLMOLCOUNT",1],
	["DBLMOLPERLITER",0,0.001],
]);

for (my $i=0; $i < @{$cdata}; $i++) {
	if ($cdata->[$i]->[3] != 0) {
		if (!defined($mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]})) {
			$mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]} = 0;
		}
		$mediacpdhash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[2]} += $cdata->[$i]->[4];
		$metahash->{$cdata->[$i]->[0]}->{$cdata->[$i]->[1]}->{associated_compounds}->{$cdata->[$i]->[2]} = $cdata->[$i]->[3]+0;
	}
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
	my $universalHash = {};
	for (my $j=0; $j < @{$universals}; $j++) {
		$universalHash->{$universals->[$j]} = 0;
	}
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
	$id = shift(@{$array});
	if (defined($array->[0]) && $array->[0] =~ /^\d+$/) {
		$id .= "_".shift(@{$array});
	}
	if (defined($idhash->{$id}) || defined($array->[0])) {
		my $index = 1;
		while (defined($idhash->{$id.".".$index})) {
			$index++;
		}
		$id = $id.".".$index;
	}
	$idhash->{$id} = 1;
	my $media = {
		id => $id,
		name => $data->[$i]->[1],
		source_id => $data->[$i]->[0],
		source => "DSMZ",
		protocol_link => $data->[$i]->[7],
		isDefined => $data->[$i]->[5]+0,
		isMinimal => 0,
		isAerobic => $data->[$i]->[6]+0,
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
		$universalHash->{$reagent} = 1;
		push(@{$media->{mediacompounds}},{
			compound_ref => "921/29/8/compounds/id/".$reagent,
			concentration => $mediacpdhash->{$data->[$i]->[0]}->{$reagent}+0,
			maxFlux => 100,
			minFlux => -100
		});
	}
	for (my $j=0; $j < @{$universals}; $j++) {
		if ($universalHash->{$universals->[$j]} == 0) {
			push(@{$media->{mediacompounds}},{
				compound_ref => "921/29/8/compounds/id/".$universals->[$j],
				concentration => 0.00001,
				maxFlux => 100,
				minFlux => -100
			});
		}
	}
	my $hasOxygen = 0;
	my $elements = {
		C => 0,
		N => 0,
		S => 0,
		P => 0
	};
	for (my $j=0; $j < @{$media->{mediacompounds}}; $j++) {
		if ($media->{mediacompounds}->[$j]->{compound_ref} =~ m/cpd00001/) {
			$media->{mediacompounds}->[$j]->{concentration} = 55;
		} elsif ($media->{mediacompounds}->[$j]->{compound_ref} =~ m/cpd00067/) {
			$media->{mediacompounds}->[$j]->{concentration} = 0.0000001;
		} elsif ($media->{mediacompounds}->[$j]->{compound_ref} =~ m/cpd00007/) {
			$hasOxygen = 1;
			if ($media->{isAerobic} == 0) {
				splice(@{$media->{mediacompounds}},$j,1);
			}
		}
		if ($media->{mediacompounds}->[$j]->{compound_ref} =~ m/(cpd\d+)/) {
			my $id = $1;
			if (defined($cpdhash->{$id}) && ($cpdhash->{$id} =~ m/C[A-Z\d]/ || $cpdhash->{$id} =~ m/C$/)) {
				$elements->{C} = 1;
			}
			if (defined($cpdhash->{$id}) && ($cpdhash->{$id} =~ m/N[A-Z\d]/ || $cpdhash->{$id} =~ m/N$/)) {
				$elements->{N} = 1;
			}
			if (defined($cpdhash->{$id}) && ($cpdhash->{$id} =~ m/S[A-Z\d]/ || $cpdhash->{$id} =~ m/S$/)) {
				$elements->{S} = 1;
			}
			if (defined($cpdhash->{$id}) && ($cpdhash->{$id} =~ m/P[A-Z\d]/ || $cpdhash->{$id} =~ m/P$/)) {
				$elements->{P} = 1;
			}
		}
	}
	if ($media->{isAerobic} == 1 && $hasOxygen == 0) {
		push(@{$media->{mediacompounds}},{
			compound_ref => "921/29/8/compounds/id/cpd00007",
			concentration => 0.000001,
			maxFlux => 20,
			minFlux => -100
		});
	}
	print $i."\t".$media->{id}."\t".$media->{name}."\t".$elements->{C}."\t".$elements->{N}."\t".$elements->{S}."\t".$elements->{P}."\n";
	#save_workspace_object("KomodoMedia/".$id,$media,"KBaseBiochem.Media");
}		