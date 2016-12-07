#!/usr/bin/perl -w

use strict;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
$|=1;

my $directory = $ARGV[0];
open ( MAIN, ">", $directory."KBaseModelList.txt");
my $list = get_ws_objects_list("ModelSEEDModels","KBaseFBA.FBAModel");
for (my $i=0; $i < @{$list}; $i++) {
	print MAIN join("\t",@{$list->[$i]})."\n";
}
close(MAIN);