#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::P3::Workspace::ScriptHelpers;
use Bio::KBase::ObjectAPI::utilities;


my $data = Bio::KBase::ObjectAPI::utilities::rest_download({url => "https://www.patricbrc.org/api/genome_feature/?genome_id=107806.33&http_accept=application/json",token => Bio::P3::Workspace::ScriptHelpers::token()});
print Data::Dumper->Dump([$data]);