#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $ws = get_ws_client();
#my $output = $ws->undelete_workspace({
#	workspace => "jimdavis:1420727110661"
#});
my $output = $ws->undelete_objects([{
	workspace => "jimdavis2:1420825381584",
	name => "Narrative.1420825381584"
}]);