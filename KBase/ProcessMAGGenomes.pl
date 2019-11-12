#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = $ARGV[0];

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $output = Bio::KBase::kbaseenv::list_objects({
	workspaces => [$workspace]
});

my $rxnhash = Bio::KBase::utilities::reaction_hash();

for (my $i=0; $i < @{$output}; $i++) {
	if ($output->[$i]->[2] =~ m/KBaseGenomes\.Genome/ && $output->[$i]->[1] =~ m/(.+)\.RAST$/) {
		system("perl ProcessMAGGenome.pl ".$output->[$i]->[1]." ".$workspace);
	}
}