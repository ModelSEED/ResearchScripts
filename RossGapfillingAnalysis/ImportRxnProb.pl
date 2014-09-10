#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(save_workspace_object printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $filename = $ARGV[0];
my $fba = get_fba_client();
my $ws = get_ws_client();

open(my $fh, "<", $filename) || return;
my $data = [];
while (my $line = <$fh>) {
	chomp($line);
	push(@{$data},[split(/\t/,$line)]);
}
close($fh);

my $genome;
my $rxnobject
for (my $i=0; $i < @{$data}; $i++) {
	if ($genome ne $data->[$i]->[0]) {
		if (defined($rxnobject)) {
			save_workspace_object("chenry:QualitativeGapfillingStudy/".$genome.".direct.rxnprob",$probobject,"ProbabilisticAnnotation.ProbAnno");
		}
		$genome = $data->[$i]->[0];
		$rxnobject = {
			id => $genome.".rxnprob",
			genome => $genome,
			genome_workspace => "chenry:QualitativeGapfillingStudy",
			template_model => "GramPosModelTemplate",
			template_workspace => "KBaseTemplateModels",
			probanno => $genome.".probanno",
			probanno_workspace => "chenry:QualitativeGapfillingStudy",
			reaction_probabilities => [],			
		};
	}
	push(@{$rxnobject->{reaction_probabilities}},[$data->[$i]->[1],$data->[$i]->[2],"NOCOMPLEXES","",""]);
}