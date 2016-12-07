#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(save_workspace_object printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_pa_client get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $filename = $ARGV[0];
my $fba = get_fba_client();
my $ws = get_ws_client();
my $pa = get_pa_client();

open(my $fh, "<", $filename) || return;
my $data = [];
while (my $line = <$fh>) {
	chomp($line);
	push(@{$data},[split(/\t/,$line)]);
}
close($fh);

my $genome;
my $probobject
for (my $i=0; $i < @{$data}; $i++) {
	if ($genome ne $data->[$i]->[0]) {
		if (defined($probobject)) {
			save_workspace_object("chenry:QualitativeGapfillingStudy/".$genome.".probanno",$probobject,"ProbabilisticAnnotation.ProbAnno");
			$pa->calculate({
				probanno => $genome.".probanno",
				probanno_workspace => "chenry:QualitativeGapfillingStudy",
				template_model => "GramPosModelTemplate",
				template_model_workspace => "KBaseTemplateModels",
				rxnprobs => $genome.".probanno.rxnprobs",
				rxnprobs_workspace => "chenry:QualitativeGapfillingStudy",
			});
		}
		$genome = $data->[$i]->[0];
		$probobject = {
			id => $genome.".probanno",
			genome => $genome,
			genome_workspace => "chenry:QualitativeGapfillingStudy",
			roleset_probabilities => {$genome.".peg.0" => []},
			skipped_features => []			
		};
	}
	push(@{$probobject->{roleset_probabilities}->{$genome.".peg.0"}},[$data->[$i]->[1],$data->[$i]->[2]]);
}