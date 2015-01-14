#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $procindex = $ARGV[0];
my $numprocs = $ARGV[1];

my $fba = get_fba_client();
my $fbaout;
my $rxnhash;
my $fbas = get_ws_objects_list("chenry:MetabliteModelingAnalysis","KBaseFBA.FBA");
for (my $i=0; $i < @{$fbas}; $i++) {
#for (my $i=0; $i < 20; $i++) {
	print $i."\n";
	(my $fbaobj,my $info) = get_workspace_object($fbas->[$i]->[6]."/".$fbas->[$i]->[0]."/".$fbas->[$i]->[4]);
	my $rxns = $fbaobj->{FBAReactionVariables};
	for (my $j=0; $j < @{$rxns}; $j++) {
		if ($rxns->[$j]->{modelreaction_ref} =~ m/(rxn\d+)/) {
			my $id = $1;
			if (abs($rxns->[$j]->{value}) > 0.0000001) {
				if (!defined($rxnhash->{$id})) {
					$rxnhash->{$id} = [0,0];
				}
				$rxnhash->{$id}->[0]++;
			}
			if ($rxns->[$j]->{max} < -0.0000001) {
				if (!defined($rxnhash->{$id})) {
					$rxnhash->{$id} = [0,0];
				}
				$rxnhash->{$id}->[1]++;
			}
			if ($rxns->[$j]->{min} > 0.0000001) {
				if (!defined($rxnhash->{$id})) {
					$rxnhash->{$id} = [0,0];
				}
				$rxnhash->{$id}->[1]++;
			}
		}
	}
}
(my $modobj,my $info) = get_workspace_object("jplfaria:modelingtranscriptomics/iBsu1103");
my $rxns = $modobj->{modelreactions};
my $cpdhash = {};
for (my $j=0; $j < @{$rxns}; $j++) {
	for (my $k=0; $k < @{$rxns->[$j]->{modelReactionReagents}}; $k++) {
		my $rxnid;
		if ($rxns->[$j]->{id} =~ m/(rxn\d+)/) {
			$rxnid = $1;
			if ($rxns->[$j]->{modelReactionReagents}->[$k]->{modelcompound_ref} =~ m/(cpd\d+)/) {
				my $id = $1;
				if (!defined($cpdhash->{$id})) {
					$cpdhash->{$id} = {
						active => 0,
						essential => 0,
						reactions => {}
					};
				}
				if (!defined($cpdhash->{$id}->{reactions}->{$rxnid})) {
					if (defined($rxnhash->{$rxnid})) {
						$cpdhash->{$id}->{active} += $rxnhash->{$rxnid}->[0];
						$cpdhash->{$id}->{essential} += $rxnhash->{$rxnid}->[1];
						if (!defined($cpdhash->{$id}->{reactions}->{$rxnid})) {
							$cpdhash->{$id}->{reactions}->{$rxnid} = 0;
						}
						$cpdhash->{$id}->{reactions}->{$rxnid}++;
					}
				}
			}
		}
	}
}
($modobj,$info) = get_workspace_object("chenrydemo/83333.1.fbamdl");
$rxns = $modobj->{modelreactions};
for (my $j=0; $j < @{$rxns}; $j++) {
	for (my $k=0; $k < @{$rxns->[$j]->{modelReactionReagents}}; $k++) {
		my $rxnid;
		if ($rxns->[$j]->{id} =~ m/(rxn\d+)/) {
			$rxnid = $1;
			if ($rxns->[$j]->{modelReactionReagents}->[$k]->{modelcompound_ref} =~ m/(cpd\d+)/) {
				my $id = $1;
				if (!defined($cpdhash->{$id})) {
					$cpdhash->{$id} = {
						active => 0,
						essential => 0,
						reactions => {}
					};
				}
				if (!defined($cpdhash->{$id}->{reactions}->{$rxnid})) {
					if (defined($rxnhash->{$rxnid})) {
						$cpdhash->{$id}->{active} += $rxnhash->{$rxnid}->[0];
						$cpdhash->{$id}->{essential} += $rxnhash->{$rxnid}->[1];
						if (!defined($cpdhash->{$id}->{reactions}->{$rxnid})) {
							$cpdhash->{$id}->{reactions}->{$rxnid} = 0;
						}
						$cpdhash->{$id}->{reactions}->{$rxnid}++;
					}
				}
			}
		}
	}
}

print "Compound\tActive count\tEssential count\tReaction count\tReactions\n";
foreach my $cpd (keys(%{$cpdhash})) {
	print $cpd."\t".$cpdhash->{$cpd}->{active}."\t".$cpdhash->{$cpd}->{essential}."\t".keys(%{$cpdhash->{$cpd}->{reactions}})."\t".join(";",keys(%{$cpdhash->{$cpd}->{reactions}}))."\n";
}	