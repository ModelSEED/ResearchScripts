#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $ess = [qw(
cpd00001
cpd00007
cpd00011
cpd00030
cpd00034
cpd00058
cpd00063
cpd00067
cpd00099
cpd00104
cpd00149
cpd00205
cpd00218
cpd00254
cpd00305
cpd00971
cpd08021
cpd10515
cpd10516
)];

my $model = $ARGV[0];
my $fba = get_fba_client();
my $medialist = get_ws_objects_list("KBaseMedia","KBaseBiochem.Media");
for (my $i=0; $i < @{$medialist}; $i++) {
	print $i."\n";
	my ($media,$info) = get_workspace_object("KBaseMedia/".$medialist->[$i]->[1]);
	my $cpds = $media->{mediacompounds};
	my $hash = {};
	my $baseref;
	for (my $j=0; $j < @{$cpds}; $j++) {
		if ($cpds->[$j]->{compound_ref} =~ m/(.+)(cpd\d+)/) {
			$baseref = $1;
			my $cpd = $2;
			if ($cpd eq "cpd00027" && $medialist->[$i]->[1] !~ m/^Carbon-/) {
				$cpd = "cpd00036";
			}
			$cpds->[$j]->{compound_ref} = $baseref.$cpd;
			$hash->{$cpd} = 1;
		}
	}
	for (my $j=0; $j < @{$ess}; $j++) {
		if (!defined($hash->{$ess->[$j]})) {
			push(@{$media->{mediacompounds}},{
				compound_ref => $baseref.$ess->[$j],
				concentration => 0.001,
				maxFlux => 100,
				minFlux => -100
			});
		}
	}
	save_workspace_object("RhodoMedia/".$medialist->[$i]->[1],$media,"KBaseBiochem.Media");
}