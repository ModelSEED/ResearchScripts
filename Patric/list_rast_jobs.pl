#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::P3::Workspace::ScriptHelpers;
use Bio::ModelSEED::MSSeedSupportServer::MSSeedSupportClient;

my $msclient = Bio::ModelSEED::MSSeedSupportServer::MSSeedSupportClient->new("http://bio-data-1.mcs.anl.gov/services/ms_fba",token => Bio::P3::Workspace::ScriptHelpers::token());
my $output = $msclient->list_rast_jobs({});
print Data::Dumper->Dump([$output]);