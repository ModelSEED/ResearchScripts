#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = $ARGV[0];
my $media = $ARGV[1];
my $template = $ARGV[2];

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $output = Bio::KBase::kbaseenv::list_objects({
	workspaces => [$workspace]
});

for (my $i=0; $i < @{$output}; $i++) {
	if ($output->[$i]->[2] =~ m/KBaseGenomes\.Genome/) {
		Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
			workspace => $workspace,
			genome_id => $output->[$i]->[1],
			fbamodel_output_id => $output->[$i]->[1].".model",
			media_id => $media,
			template_id => $template,
			genome_workspace => $workspace,
			template_workspace => $workspace,
			media_workspace => $workspace,
			gapfill_model => 1,
		});
	}
}