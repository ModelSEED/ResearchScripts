#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $procindex = $ARGV[0];
my $numprocs = $ARGV[1];
my $fba = get_fba_client();

my $genomes = get_ws_objects_list("KBasePublicGenomesV4","KBaseGenomes.Genome");

for (my $i=0; $i < @{$genomes}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		my $output = $fba->genome_to_fbamodel({
			genome_workspace => "KBasePublicGenomesV4",
			genome => $genomes->[$i]->[1],
			workspace => "KBasePublicModelsV4",
			model => $genomes->[$i]->[1].".model"
	   	});
		printObjectInfo($output);
		$output = $fba->gapfill_model({
			model => $genomes->[$i]->[1].".model",
			workspace => "KBasePublicModelsV4",
			out_model => $genomes->[$i]->[1].".model.gf",
			solver => "cplex",
			fastgapfill => 1,
			integrate_solution => 1
		});
		printObjectInfo($output);
	}
}


		