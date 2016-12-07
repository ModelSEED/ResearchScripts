#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::P3::Workspace::ScriptHelpers;

my $usage = "Usage:\nget_node <Node ID> <filename> <Service URL>\n";
if ((defined($ARGV[0]) && $ARGV[0] eq "-h") || !defined($ARGV[0])) {
	print $usage;
	exit 0;
}
if (!defined($ARGV[2])) {
	$ARGV[2] = "http://p3.theseed.org/services/shock_api";
}

system('curl -X GET -H "Authorization: OAuth '.Bio::P3::Workspace::ScriptHelpers::token().'" '.$ARGV[2].'/node/'.$ARGV[0].'?download > '.$ARGV[1]);