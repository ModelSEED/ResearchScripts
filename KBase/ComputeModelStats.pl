#!/usr/bin/perl -w

use strict;
use Bio::KBase::workspace::ScriptHelpers qw(workspaceURL get_ws_client workspace parseObjectMeta parseWorkspaceMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(getToken fbaws get_fba_client runFBACommand universalFBAScriptCode );

$|=1;

my $fba = get_fba_client();
my $ws = get_ws_client();
my $wsmodels = $ws->list_objects({
	workspaces => ["HL_PNNLGeneCallsKBaseAnnotation"],
	type => "KBaseFBA.FBAModel"
});

my $top = [qw(
		id
		scientific_name
		dna_size
		num_contigs
		genome_ref
		domain
		taxonomy
		source
		gc_content
		total_reactions
		total_genes
		total_compounds
		extracellular_compounds
		intracellular_compounds
		transport_reactions
		subsystem_reactions
		subsystem_genes
		spontaneous_reactions
		reactions_with_genes
		gapfilled_reactions
 		model_genes
		minimal_essential_genes
		complete_essential_genes
		minimal_essential_reactions
		complete_essential_reactions
		minimal_blocked_reactions
		complete_blocked_reactions
		minimal_variable_reactions
		complete_variable_reactions
		growth_complete_media
		growth_minimal_media
)];

my $bottom = [qw(
		class
		subclass
		genes
		reactions
		model_genes
		minimal_essential_genes
		complete_essential_genes
		minimal_essential_reactions
		complete_essential_reactions
		minimal_blocked_reactions
		complete_blocked_reactions
		minimal_variable_reactions
		complete_variable_reactions
)];

my $data;
my $ssdata;
my $models;
for (my $i=0; $i < @{$wsmodels}; $i++) {
	if ($wsmodels->[$i]->[1] =~ m/model_gapfilled_transporters_added/) {
		print STDERR "Processing model:".$wsmodels->[$i]->[1]."\n";
		my $output = $fba->generate_model_stats({
			model => $wsmodels->[$i]->[1],
			workspace => "HL_PNNLGeneCallsKBaseAnnotation"
		});
		push(@{$models},$wsmodels->[$i]->[1]);
		for (my $j=0; $j < @{$top}; $j++) {
			$data->{$top->[$j]}->{$wsmodels->[$i]->[1]} = $output->{$top->[$j]};
		}
		foreach my $ss (@{$output->{subsystems}}) {
			for (my $j=0; $j < @{$bottom}; $j++) {
				$ssdata->{$ss->{name}}->{$bottom->[$j]}->{$wsmodels->[$i]->[1]} = $ss->{$bottom->[$j]};
			}
		}
	}
}
print STDERR "Processing model:bin16_OSCAR_PcalledPannotated_model_gapfilled_transporters_added\n";
my $output = $fba->generate_model_stats({
	model => "bin16_OSCAR_PcalledPannotated_model_gapfilled_transporters_added",
	workspace => "HL_PNNLGeneCallsAndAnnotation"
});
push(@{$models},"bin16_OSCAR_PcalledPannotated_model_gapfilled_transporters_added");
for (my $j=0; $j < @{$top}; $j++) {
	$data->{$top->[$j]}->{"bin16_OSCAR_PcalledPannotated_model_gapfilled_transporters_added"} = $output->{$top->[$j]};
}
foreach my $ss (@{$output->{subsystems}}) {
	for (my $j=0; $j < @{$bottom}; $j++) {
		$ssdata->{$ss->{name}}->{$bottom->[$j]}->{"bin16_OSCAR_PcalledPannotated_model_gapfilled_transporters_added"} = $ss->{$bottom->[$j]};
	}
}
print STDERR "Processing model:bin11_ANA_PcalledPannotated_model_gapfilled_transporters_added\n";
$output = $fba->generate_model_stats({
	model => "bin11_ANA_PcalledPannotated_model_gapfilled_transporters_added",
	workspace => "HL_PNNLGeneCallsAndAnnotation"
});
push(@{$models},"bin11_ANA_PcalledPannotated_model_gapfilled_transporters_added");
for (my $j=0; $j < @{$top}; $j++) {
	$data->{$top->[$j]}->{"bin11_ANA_PcalledPannotated_model_gapfilled_transporters_added"} = $output->{$top->[$j]};
}
foreach my $ss (@{$output->{subsystems}}) {
	for (my $j=0; $j < @{$bottom}; $j++) {
		$ssdata->{$ss->{name}}->{$bottom->[$j]}->{"bin11_ANA_PcalledPannotated_model_gapfilled_transporters_added"} = $ss->{$bottom->[$j]};
	}
}

print "Printing model summary data:\n";
print "Models\t".join("\t",@{$models})."\n";
for (my $j=0; $j < @{$top}; $j++) {
	print $top->[$j];
	for (my $i=0; $i < @{$models}; $i++) {
		if (defined($data->{$top->[$j]}->{$models->[$i]})) {
			print "\t".$data->{$top->[$j]}->{$models->[$i]};
		} else {
			print "\t";
		}
	}
	print "\n";
}
print "\nPrinting SS summary data:\n";
for (my $i=0; $i < @{$models}; $i++) {
	print "Subsystems";
	for (my $j=0; $j < @{$bottom}; $j++) {
		print "\t".$models->[$i]."-".$bottom->[$j];
	}
}
print "\n";
foreach my $ss (keys(%{$ssdata})) {
	print $ss;
	for (my $i=0; $i < @{$models}; $i++) {
		for (my $j=0; $j < @{$bottom}; $j++) {
			if (defined($ssdata->{$ss}->{$bottom->[$j]}->{$models->[$i]})) {
				print "\t".$ssdata->{$ss}->{$bottom->[$j]}->{$models->[$i]};
			} else {
				print "\t";
			}
		}
	}
	print "\n";
}