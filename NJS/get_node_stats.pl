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

my $tokenObj = Bio::KBase::AuthToken->new(
    user_id => "reviewer", password => 'reviewer',ignore_authrc => 1
);
my $token = $tokenObj->token();
my $node = "73ac5fa8-cc41-4713-9678-0b7b6a147ae8";
my $shockurl = "http://p3.theseed.org/services/shock_api/node/".$node;
my $ua = LWP::UserAgent->new();
my $res = $ua->get($shockurl,Authorization => "OAuth ".$token);
#print Data::Dumper->Dump([$res]);
my $json = JSON::XS->new;
my $data = $json->decode($res->content);
if (defined($data->{data}->{file}->{size})) {
	print "Size = ".$data->{data}->{file}->{size}."\n";
}
#print "create shock node output:\n".Data::Dumper->Dump([$data])."\n\n";


