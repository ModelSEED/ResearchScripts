#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(getToken fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $usage = "Usage:\nget_node <Node ID> <Service URL> <filename>\n";
if ((defined($ARGV[0]) && $ARGV[0] eq "-h") || !defined($ARGV[0])) {
	print $usage;
	exit 0;
}
if (!defined($ARGV[1])) {
	$ARGV[1] = "http://kbase.us/services/shock-api";
}

system('curl -X GET -H "Authorization: OAuth '.getToken().'" '.$ARGV[1].'/node/'.$ARGV[0].'?download > '.$ARGV[2]);