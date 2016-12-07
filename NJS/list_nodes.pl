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

my $usage = "Usage:\nlist_nodes <Service URL>\n";
if ((defined($ARGV[0]) && $ARGV[0] eq "-h")) {
	print $usage;
	exit 0;
}
if (!defined($ARGV[0])) {
	$ARGV[0] = "http://kbase.us/services/shock-api";
}

my $var = 'curl -X GET '.$ARGV[0]."/node?limit=25000";
system($var);