#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $damagergts = [qw(
cpd00229
cpd00428
cpd01160
cpd00215
cpd00146
cpd00084
cpd00547
cpd00587
cpd10147
cpd03604
cpd17158
cpd01809
cpd01495
cpd04254
cpd12138
cpd00216
cpd00658
cpd01017
cpd15082
cpd00040
cpd05186
cpd11311
cpd00219
cpd02951
cpd02980
cpd00616
cpd00231
cpd00493
cpd02978
cpd00095
cpd01671
cpd00157
cpd04532
cpd00016
cpd06805
cpd03603
cpd17397
cpd03530
cpd03666
cpd00858
cpd00434
cpd12219
cpd00340
cpd00053
cpd00253
cpd00610
cpd00814
cpd02187
cpd03484
cpd00199
cpd00056
cpd00305
cpd01852
cpd00793
cpd01937
cpd11180
cpd00004
cpd00005
cpd01296
cpd00015
cpd00982
cpd00584
cpd01775
cpd01777
cpd02820
cpd01140
cpd02177
cpd02338
cpd00103
cpd00932
cpd16335
cpd01418
cpd16381
cpd12176
cpd03518
cpd11231
cpd00957
cpd03519
cpd14621
cpd14622
cpd14730
cpd14859
cpd00931
cpd02642
cpd02508
cpd02728
cpd03521
cpd00352
cpd02431
cpd02651
cpd01793
cpd00236
cpd00817
cpd01994
cpd02081
cpd04185
cpd17522
cpd00896
cpd03606
cpd01895
cpd02961
cpd02186
cpd01112
cpd01152
cpd15380
cpd15494
cpd16540
cpd01875
cpd02122
cpd03889
cpd11314
cpd08636
cpd00551
cpd00449
cpd02315
cpd01454
cpd01872
cpd07807
cpd00884
cpd02394
cpd17953
cpd03470
cpd14545
cpd14905
cpd16322
cpd16376
cpd03411
cpd16380
cpd02443
cpd02678
cpd02826
cpd00706
cpd00059
cpd00766
)];
my $damagehash = {};
for (my $i=0; $i < @{$damagergts}; $i++) {
	$damagehash->{$damagergts->[$i]} = 1;
}
my $rxnhash;
my $mdlhash;
my $mdlobjhash;
my $fbas = get_ws_objects_list("MetaboliteEssentialityAnalysis","KBaseFBA.FBA");
for (my $i=0; $i < @{$fbas}; $i++) {
	print $i."\n";
	(my $fbaobj,my $info) = get_workspace_object($fbas->[$i]->[6]."/".$fbas->[$i]->[0]."/".$fbas->[$i]->[4]);
	my $ref = $fbaobj->{fbamodel_ref};
	if (!defined($mdlobjhash->{$ref})) {
		($mdlobjhash->{$ref},$info) = get_workspace_object($ref);
		my $modobj = $mdlobjhash->{$ref};
		my $rxns = $modobj->{modelreactions};
		for (my $j=0; $j < @{$rxns}; $j++) {
			my $rxnid = $rxns->[$j]->{id};
			$mdlobjhash->{$ref}->{rxnhash}->{$rxnid} = $rxns->[$j];
			for (my $k=0; $k < @{$rxns->[$j]->{modelReactionReagents}}; $k++) {
				if ($rxns->[$j]->{modelReactionReagents}->[$k]->{modelcompound_ref} =~ m/\/([^\/]+)$/) {
					my $id = $1;
					$rxns->[$j]->{rgts}->{$id} = $rxns->[$j]->{modelReactionReagents}->[$k]->{coefficient};
					$modobj->{cpdhash}->{$id}->{$rxnid} = $rxns->[$j]->{modelReactionReagents}->[$k]->{coefficient};
					my $array = [split(/_/,$id)];
					$id = $array->[0];
					if (defined($damagehash->{$id})) {
						$rxns->[$j]->{damagecpd}->{$id} = 1;
					}
				}
			}
		}
	}
	my $rxns = $fbaobj->{FBAReactionVariables};
	for (my $j=0; $j < @{$rxns}; $j++) {
		if ($rxns->[$j]->{modelreaction_ref} =~ m/\/([^\/]+)$/) {
			$fbaobj->{rxnhash}->{$1} = $rxns->[$j];
		}
	}
	for (my $j=0; $j < @{$rxns}; $j++) {
		if ($rxns->[$j]->{modelreaction_ref} =~ m/\/([^\/]+)$/) {
			my $id = $1;
			my $dir = "F";
			my $rxnobj = $mdlobjhash->{$ref}->{rxnhash}->{$id};
			my $dmgrgt = keys(%{$rxnobj->{damagecpd}});
			my $sidergt = 0;
			my $sideessrgt = 0;
			foreach my $rgt (keys(%{$rxnobj->{rgts}})) {
				my $rxncount = keys(%{$mdlobjhash->{cpdhash}->{$rgt}});
				if ($rxncount < 15) {
					my $found = 0;
					my $foundess = 0;
					foreach my $rxn (keys(%{$mdlobjhash->{cpdhash}->{$rgt}})) {
						if ($rxn ne $id && defined($fbaobj->{rxnhash}->{$rxn})) {
							if ($rxns->[$j]->{max} < -0.0000001 || $rxns->[$j]->{min} > 0.0000001) {
								$foundess = 1;
							} elsif ($fbaobj->{rxnhash}->{$rxn}->{max} > 0.0000001 || $fbaobj->{rxnhash}->{$rxn}->{min} < 0.0000001) {
								$found = 1;
							}
						}
					}
					if ($found == 1) {
						$sidergt++;
					}
					if ($foundess == 1) {
						$sideessrgt++;
					}
				}
			}
			if (abs($rxns->[$j]->{value}) > 0.0000001) {
				if ($rxns->[$j]->{value} < 0) {
					$dir = "R";
				}
				if (!defined($rxnhash->{$id}->{$dir})) {
					$rxnhash->{$id}->{$dir} = {
						active => 0,
						essential => 0,
						activemodels => 0,
						essentialmodels => 0,
						aveessflux => 0,
						aveactflux => 0,
						damagereactants => $dmgrgt,
						sidereactants => $sidergt,
						esssides => $sideessrgt
					};
				}
				if (!defined($mdlhash->{$id}->{$dir}->{$ref})) {
					$mdlhash->{$id}->{$dir}->{$ref} = 1;
					$rxnhash->{$id}->{$dir}->{activemodels}++;
				}
				$rxnhash->{$id}->{$dir}->{active}++;
				$rxnhash->{$id}->{$dir}->{aveactflux} += abs($rxns->[$j]->{value});
			}
			if ($rxns->[$j]->{max} < -0.0000001 || $rxns->[$j]->{min} > 0.0000001) {
				if (!defined($rxnhash->{$id}->{$dir})) {
					$rxnhash->{$id}->{$dir} = {
						active => 0,
						essential => 0,
						activemodels => 0,
						essentialmodels => 0,
						aveessflux => 0,
						aveactflux => 0,
						damagereactants => $dmgrgt,
						sidereactants => $sidergt,
						esssides => $sideessrgt
					};
				}
				if (!defined($mdlhash->{$id}->{$dir}->{$ref})) {
					$mdlhash->{$id}->{$dir}->{$ref} = 1;
					$rxnhash->{$id}->{$dir}->{essentialmodels}++;
				}
				$rxnhash->{$id}->{$dir}->{essential}++;
				$rxnhash->{$id}->{$dir}->{aveessflux} += abs($rxns->[$j]->{value});
			}
		}
	}
}
foreach my $key (keys(%{$rxnhash})) {
	foreach my $dir (keys(%{$rxnhash->{$key}})) {
		if ($rxnhash->{$key}->{$dir}->{essential} > 0) {
			$rxnhash->{$key}->{$dir}->{aveessflux} = $rxnhash->{$key}->{$dir}->{aveessflux}/$rxnhash->{$key}->{$dir}->{essential};
		}
		if ($rxnhash->{$key}->{$dir}->{active} > 0) {
			$rxnhash->{$key}->{$dir}->{aveessflux} = $rxnhash->{$key}->{$dir}->{aveessflux}/$rxnhash->{$key}->{$dir}->{active};
		}
	}
}
open (my $out, ">", $directory."/ReactionActivityAnalysis.txt");
print $out "Reaction\tDirection\tActive\tEssential\tActive models\tEssential models\tReactants with active side pathways\tReactants with essential side pathways\tAverage essential flux\tAverage active flux\tDamage prone reactants\n";
foreach my $key (keys(%{$rxnhash})) {
	foreach my $dir (keys(%{$rxnhash->{$key}})) {
		print $out $key."\t".$dir."\t".$rxnhash->{$key}->{$dir}->{active}."\t".$rxnhash->{$key}->{$dir}->{essential}."\t".$rxnhash->{$key}->{$dir}->{activemodels}."\t".$rxnhash->{$key}->{$dir}->{essentialmodels}."\t".$rxnhash->{$key}->{$dir}->{sidereactants}."\t".$rxnhash->{$key}->{$dir}->{esssides}."\t".$rxnhash->{$key}->{$dir}->{aveessflux}."\t".$rxnhash->{$key}->{$dir}->{aveactflux}."\t".$rxnhash->{$key}->{$dir}->{damagereactants}."\n"
	}
}
close($out);