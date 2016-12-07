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

my $ws = "chenry:1440386682869";
my $rxnhash = {};
my $cpdhash = {};
my $genehash = {};
(my $data,my $prov) = get_workspace_object($ws."/SB2B_model.gf");
my $mdlrxns = $data->{modelreactions};
for (my $i=0; $i < @{$mdlrxns}; $i++) {
	my $id;
	if ($mdlrxns->[$i]->{reaction_ref} =~ m/(rxn\d+)$/) {
		if ($1 ne "rxn00000") {
			$id = $1;
		}
	}
	if (defined($id) && $mdlrxns->[$i]->{modelcompartment_ref} =~ m/([a-z]\d+)$/) {
		$id .= "_".$1;
	}
	if (!defined($id)) {
		$id = $mdlrxns->[$i]->{id};
	}
	$rxnhash->{$id} = [$i,0];
	my $prots = $mdlrxns->[$i]->{modelReactionProteins};
	for (my $j=0; $j < @{$prots}; $j++) {
		my $sus = $prots->[$j]->{modelReactionProteinSubunits};
		for (my $k=0; $k < @{$sus}; $k++) {
			my $ftrs = $sus->[$k]->{feature_refs};
			for (my $m=0; $m < @{$ftrs}; $m++) {
				if ($ftrs->[$m] =~ m/\/([^\/]+$)/) {
					my $gene = $1;
					$genehash->{$gene} = 1;
					$rxnhash->{$id} = [$i,1];
				}
			}
		}
	}
}
my $mdlcpds = $data->{modelcompounds};
for (my $i=0; $i < @{$mdlcpds}; $i++) {
	$cpdhash->{$mdlcpds->[$i]->{id}} = $mdlcpds->[$i];
}

(my $datat,$prov) = get_workspace_object($ws."/SB2B_draft_model");
my $othercpdhash = {};
$mdlcpds = $datat->{modelcompounds};
for (my $i=0; $i < @{$mdlcpds}; $i++) {
	$othercpdhash->{$mdlcpds->[$i]->{id}} = $mdlcpds->[$i];
}
$mdlrxns = $datat->{modelreactions};
for (my $i=0; $i < @{$mdlrxns}; $i++) {
	if (defined($rxnhash->{$mdlrxns->[$i]->{id}}) && $rxnhash->{$mdlrxns->[$i]->{id}}->[1] == 0) {
		$data->{modelreactions}->[$rxnhash->{$mdlrxns->[$i]->{id}}->[0]]->{modelReactionProteins} = $mdlrxns->[$i]->{modelReactionProteins};
	} elsif (!defined($rxnhash->{$mdlrxns->[$i]->{id}})) {
		my $add = 0;
		my $prots = $mdlrxns->[$i]->{modelReactionProteins};
		for (my $j=0; $j < @{$prots}; $j++) {
			my $sus = $prots->[$j]->{modelReactionProteinSubunits};
			for (my $k=0; $k < @{$sus}; $k++) {
				my $ftrs = $sus->[$k]->{feature_refs};
				for (my $m=0; $m < @{$ftrs}; $m++) {
					if ($ftrs->[$m] =~ m/\/([^\/]+$)/) {
						my $gene = $1;
						if ($gene =~ m/SAMA/ && !defined($genehash->{$gene})) {
							$add = 1;	
						}
					}
				}
			}
		}
		if ($add == 1) {
			push(@{$data->{modelreactions}},$mdlrxns->[$i]);
			my $rgts = $mdlrxns->[$i]->{modelReactionReagents};
			for (my $j=0; $j < @{$rgts}; $j++) {
				if ($rgts->[$j]->{modelcompound_ref} =~ m/(cpd.+)$/) {
					my $id = $1;
					if (!defined($cpdhash->{$id})) {
						$cpdhash->{$id} = $othercpdhash->{$id};
						push(@{$data->{modelcompounds}},$cpdhash->{$id});
					}
				}
			}
		}
	}
}
save_workspace_object($ws."/SB2B_reconciled_model",$data,"KBaseFBA.FBAModel");