#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];
my $module = $ARGV[1];

open (my $fh, "<", "/Users/chenry/workspace/KBasePythonModules/".$module);
my $code = "";
while (my $line = <$fh>) {
    $code .= $line;
}
close($fh);

my $obj = {
	description => $code
};

save_workspace_object("KBasePublicModules/".$module,$obj,"KBaseNarrative.Metadata");
