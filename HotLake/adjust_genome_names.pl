#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $ws = fbaws();

my $genomes = [qw(
bin17-genome-PcalledKannotated
bin10-genome-PcalledKannotated
bin19-genome-PcalledKannotated
bin12-genome-PcalledKannotated
bin05-genome-PcalledKannotated
bin07-genome-PcalledKannotated
bin14-genome-PcalledKannotated
bin13-genome-PcalledKannotated
bin24-genome-PcalledKannotated
bin02-genome-PcalledKannotated
bin06-genome-PcalledKannotated
bin21-genome-PcalledKannotated
bin20-genome-PcalledKannotated
bin03-genome-PcalledKannotated
bin18-genome-PcalledKannotated
bin15-genome-PcalledKannotated
bin08-genome-PcalledKannotated
bin09-genome-PcalledKannotated
bin04-genome-PcalledKannotated
bin01-genome-PcalledKannotated
bin22-genome-PcalledKannotated
)];

for (my $i=0; $i < @{$genomes}; $i++) {
	(my $data,my $prov) = get_workspace_object($ws."/".$genomes->[$i]);
	if ($genomes->[$i] =~ m/(bin\d+)/) {
		$data->{scientific_name} = $1." genome";
	}
	save_workspace_object($ws."/".$genomes->[$i],$data,"KBaseGenomes.Genome");
}
