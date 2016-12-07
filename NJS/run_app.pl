#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use JSON::XS;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(getToken fbaws get_fba_client runFBACommand universalFBAScriptCode );
use Bio::KBase::NarrativeJobService::Client;
my $usage = "Usage:\nTestNJS <Input file> <Service URL>\n";
if ((defined($ARGV[0]) && $ARGV[0] eq "-h") || !defined($ARGV[0])) {
	print $usage;
	exit 0;
}
if (!defined($ARGV[1])) {
	$ARGV[1] = "http://140.221.66.246:7080";
}
my $njs_obj = new Bio::KBase::NarrativeJobService::Client($ARGV[1],token => getToken());
#Loading parameters from file
open( my $fh, "<", $ARGV[0]);
my $app;
{
    local $/;
    my $str = <$fh>;
    $app = decode_json $str;
}
close($fh);
my $appstate = $njs_obj->run_app($app);
print "Output:".Data::Dumper->Dump([$appstate])."\n";