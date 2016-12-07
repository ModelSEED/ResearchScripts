#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $fba = get_fba_client();
my $fbaout;
my $rxnhash;
my $fbas = get_ws_objects_list("MetaboliteEssentialityAnalysis","KBaseFBA.FBA");
print @{$fbas}." fbas\n";
for (my $i=0; $i < @{$fbas}; $i++) {
	print $i."\n";
	(my $fbaobj,my $info) = get_workspace_object($fbas->[$i]->[6]."/".$fbas->[$i]->[0]."/".$fbas->[$i]->[4]);
	my $ref = $fbaobj->{fbamodel_ref};
	my $rxns = $fbaobj->{FBAReactionVariables};
	for (my $j=0; $j < @{$rxns}; $j++) {
		if ($rxns->[$j]->{modelreaction_ref} =~ m/\/([^\/]+)$/) {
			my $id = $1;
			if (abs($rxns->[$j]->{value}) > 0.0000001) {
				if (!defined($rxnhash->{$ref}->{$id})) {
					$rxnhash->{$ref}->{$id} = [0,0];
				}
				$rxnhash->{$ref}->{$id}->[0]++;
			}
			if ($rxns->[$j]->{max} < -0.0000001) {
				if (!defined($rxnhash->{$ref}->{$id})) {
					$rxnhash->{$ref}->{$id} = [0,0];
				}
				$rxnhash->{$ref}->{$id}->[1]++;
			}
			if ($rxns->[$j]->{min} > 0.0000001) {
				if (!defined($rxnhash->{$ref}->{$id})) {
					$rxnhash->{$ref}->{$id} = [0,0];
				}
				$rxnhash->{$ref}->{$id}->[1]++;
			}
		}
	}
}
my $cpdhash = {};
foreach my $key (keys(%{$rxnhash})) {
	(my $modobj,my $info) = get_workspace_object($key);
	my $rxns = $modobj->{modelreactions};
	for (my $j=0; $j < @{$rxns}; $j++) {
		for (my $k=0; $k < @{$rxns->[$j]->{modelReactionReagents}}; $k++) {
			my $rxnid = $rxns->[$j]->{id};
			if ($rxns->[$j]->{modelReactionReagents}->[$k]->{modelcompound_ref} =~ m/\/([^\/]+)$/) {
				my $id = $1;
				if (!defined($cpdhash->{$id})) {
					$cpdhash->{$id} = {
						reactions => {}
					};
				}
				if (!defined($cpdhash->{$id}->{$key})) {
					$cpdhash->{$id}->{$key} = {
						active => 0,
						essential => 0,
						reactions => {}
					};
				}
				if (!defined($cpdhash->{$id}->{$key}->{reactions}->{$rxnid})) {
					if (defined($rxnhash->{$key}->{$rxnid})) {
						$cpdhash->{$id}->{$key}->{active} += $rxnhash->{$key}->{$rxnid}->[0];
						$cpdhash->{$id}->{$key}->{essential} += $rxnhash->{$key}->{$rxnid}->[1];
						if (!defined($cpdhash->{$id}->{$key}->{reactions}->{$rxnid})) {
							$cpdhash->{$id}->{$key}->{reactions}->{$rxnid} = 0;
						}
						$cpdhash->{$id}->{$key}->{reactions}->{$rxnid}++;
						$cpdhash->{$id}->{reactions}->{$rxnid}++;
					}
				}
			}
		}
	}
}

print "Compound\tReactions";
foreach my $key (keys(%{$rxnhash})) {
	print "\tActive ".$key."\tEssential ".$key."\tReactions ".$key;
}
print "\n";
foreach my $cpd (keys(%{$cpdhash})) {
	print $cpd."\t".keys(%{$cpdhash->{$cpd}->{reactions}});
	foreach my $key (keys(%{$rxnhash})) {
		if (!defined($cpdhash->{$cpd}->{$key})) {
			$cpdhash->{$cpd}->{$key} = {
				active => 0,
				essential => 0,
				reactions => {}
			};
		}
		print "\t".$cpdhash->{$cpd}->{$key}->{active}."\t".$cpdhash->{$cpd}->{$key}->{essential}."\t".keys(%{$cpdhash->{$cpd}->{$key}->{reactions}});
	}
	print "\n";
}	