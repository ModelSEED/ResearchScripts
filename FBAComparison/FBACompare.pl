#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];

my $fbas = get_ws_objects_list($workspace,"KBaseFBA.FBA");
my $hash = {};
print "Reactions";
for (my $i=0; $i < @{$fbas}; $i++) {
	print "\t".$fbas->[$i]->[1];
	(my $obj,my $info) = get_workspace_object($workspace."/".$fbas->[$i]->[1]);
	my $rxns = $obj->{FBAReactionVariables};
	for (my $j=0; $j < @{$rxns}; $j++) {
		$hash->{$rxns->[$j]->{modelreaction_ref}}->[$i] = $rxns->[$j]->{value};
	}
}
print "\n";

foreach my $rxn (keys(%{$hash})) {
	print $rxn."\t".join("\t",@{$hash->{$rxn}})."\n";
}