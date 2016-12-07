#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::P3::Workspace::ScriptHelpers;

system('curl -X GET -H "Authorization: OAuth '.Bio::P3::Workspace::ScriptHelpers::token().'" '.$ARGV[0].'?download > '.$ARGV[1]);