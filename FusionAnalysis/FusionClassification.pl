#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(print_file load_file get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
use Bio::KBase::CDMI::CDMIClient;
use Bio::KBase::GenomeAnnotation::Client;

my $fba = get_fba_client();
my $cdmi = Bio::KBase::CDMI::CDMIClient->new_for_script();
my $sapsvr = ModelSEED::Client::SAP->new();
my $directory = $ARGV[0];
my $data = load_file($directory."/AllFusions.txt");

#Getting hashes of roles to reactions and subsystems
$fba->_setContext(undef,{},undef,undef,undef);
my $template = $fba->_get_msobject("ModelTemplate","KBaseTemplateModels","GramNegModelTemplate");
my $RoleToRxn = $template->simple_role_reaction_hash();
my $RoleToSubsys = $template->mapping()->RoleNameSubsystemHash();
foreach my $key (keys(%{$RoleToRxn})) {
	my $searchrole = Bio::KBase::ObjectAPI::utilities::convertRoleToSearchRole($key);
	$RoleToRxn->{$searchrole} = $RoleToRxn->{$key};
}
foreach my $key (keys(%{$RoleToSubsys})) {
	my $searchrole = Bio::KBase::ObjectAPI::utilities::convertRoleToSearchRole($key);
	$RoleToSubsys->{$searchrole} = $RoleToSubsys->{$key};
}
#Getting reactions neighbors
my $rxnneighbors = $template->biochemistry()->neighboring_reaction_hash();
my $eqns;
my $rxns = $template->biochemistry()->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$eqns->{$rxns->[$i]->id()} = $rxns->[$i]->definition();
}
my $cpdnames;
my $cpds = $template->biochemistry()->compounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$cpdnames->{$cpds->[$i]->id()} = $cpds->[$i]->name();
}

my $genes;
my $split;
my $GeneFusionData;
my $UniqueRoleCombinations;
my $NMM = [];
my $NNMM = [];
my $SMSS = [];
my $SMNS = [];
my $NMTS = [];
my $NMSS = [];
my $NMNS = [];
for (my $i=0; $i < @{$data}; $i++) {
	my $item = $data->[$i];
	my $array = [split(/\t/,$item)];
	push(@{$genes},$array->[0]);
	push(@{$split},$array->[1]);
	if (@{$genes} >= 1000) {
  		print $i." of ".@{$data}."\n";
  		my $fidDataHash = $cdmi->fids_to_feature_data($genes);
  		my $output = $cdmi->get_entity_Feature(
			$genes,
			["id","source_id"]
		);
		#Retrieving SEED functions
		my $seedgenes = [];
		for (my $j=0; $j < @{$genes}; $j++) {
			push(@{$seedgenes},$output->{$genes->[$j]}->{"source_id"});
		}
		my $functions = $sapsvr->ids_to_functions({-ids => $seedgenes});
		#Analyzing and classifying genes
  		for (my $j=0; $j < @{$genes}; $j++) {
  			my $seedid = $seedgenes->[$j];
  			$fidDataHash->{$genes->[$j]}->{feature_function} =~ s/\s*#.+//;
  			my $fnarray = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$fidDataHash->{$genes->[$j]}->{feature_function}))];
  			my $combo = join("\t",@{$fnarray});
  			my $index;
  			if (defined($UniqueRoleCombinations->{$combo})) {
  				$UniqueRoleCombinations->{$combo}->{Count}++;
  				$index = @{$UniqueRoleCombinations->{$combo}->{SEEDGenes}};
  				push(@{$UniqueRoleCombinations->{$combo}->{SEEDGenes}},$seedid);
  				push(@{$UniqueRoleCombinations->{$combo}->{KBGenes}},$genes->[$j]);
  				push(@{$UniqueRoleCombinations->{$combo}->{Sources}},"kb");
  			} else {
  				$index = 0;
  				$UniqueRoleCombinations->{$combo} = {
  					Count => 1,
  					KBRep => $genes->[$j],
					SEEDRep => $seedid,
  					KBGenes => [$genes->[$j]],
					SEEDGenes => [$seedid],
					Sources => ["kb"],
					Shared => {},
					Reactions => [],
					Roles => $fnarray,
					SS => []
  				};
  				$UniqueRoleCombinations->{$combo}->{RoleCount} = @{$fnarray};
  				my $MetabolicCount = 0;
  				my $SubsysCount = 0;
  				my $UnannotatedCount = 0;
  				my $Neighbors = 0;
  				for (my $k=0; $k < @{$fnarray}; $k++) {
  					$UniqueRoleCombinations->{$combo}->{SS}->[$k] = [];
  					$UniqueRoleCombinations->{$combo}->{Reactions}->[$k] = [];
  					my $searchname = Bio::KBase::ObjectAPI::utilities::convertRoleToSearchRole($fnarray->[$k]);
  					if (defined($RoleToSubsys->{$searchname})) {
  						$UniqueRoleCombinations->{$combo}->{SS}->[$k] = [keys(%{$RoleToSubsys->{$searchname}})];
  					}
  					if (defined($RoleToRxn->{$searchname})) {
  						my $rxns = [keys(%{$RoleToRxn->{$searchname}})];
  						for (my $m=0; $m < @{$rxns}; $m++) {
  							$UniqueRoleCombinations->{$combo}->{Reactions}->[$k]->[$m] = $rxns->[$m].":".$eqns->{$rxns->[$m]};
  						}
  						$MetabolicCount++;
  					} elsif (defined($RoleToSubsys->{$searchname})) {
  						$SubsysCount++;
  					} else {
  						$UnannotatedCount++;
  					}
  				}
  				for (my $k=0; $k < @{$fnarray}; $k++) {
  					for (my $m=$k+1; $m < @{$fnarray}; $m++) {
  						for (my $n=0; $n < @{$UniqueRoleCombinations->{$combo}->{Reactions}->[$k]}; $n++) {
  							my $rxnone = $UniqueRoleCombinations->{$combo}->{Reactions}->[$k]->[$n];
  							$rxnone =~ s/:.+//;
  							for (my $o=0; $o < @{$UniqueRoleCombinations->{$combo}->{Reactions}->[$m]}; $o++) {
  								my $rxntwo = $UniqueRoleCombinations->{$combo}->{Reactions}->[$m]->[$o];
  								$rxntwo =~ s/:.+//;
  								if (defined($rxnneighbors->{$rxnone}->{$rxntwo}) && $rxnone ne $rxntwo) {
  									foreach my $cpd (keys(%{$rxnneighbors->{$rxnone}->{$rxntwo}})) {
 										$Neighbors++;
 										my $cpdarr = [split(/_/,$cpd)];
  										$UniqueRoleCombinations->{$combo}->{Shared}->{$cpdnames->{$cpdarr->[0]}} = 1;	
  									}
  								}
  							}
  						}
  					}
  				}
  				if ($MetabolicCount > 0) {
  					if ($MetabolicCount > 1) {
  						if ($Neighbors > 0) {
  							push(@{$NMM},$UniqueRoleCombinations->{$combo});
  						} else {
  							push(@{$NNMM},$UniqueRoleCombinations->{$combo});
  						}
  					} elsif ($SubsysCount > 0) {
  						push(@{$SMSS},$UniqueRoleCombinations->{$combo});
  					} else {
  						push(@{$SMNS},$UniqueRoleCombinations->{$combo});
  					}
  				} elsif ($SubsysCount > 0) {
  					if ($SubsysCount > 1) {
  						push(@{$NMTS},$UniqueRoleCombinations->{$combo});
  					} else {
  						push(@{$NMSS},$UniqueRoleCombinations->{$combo});
  					}
  				} else {
  					push(@{$NMNS},$UniqueRoleCombinations->{$combo});
  				}
  			}
  			if (defined($seedid)) {
	  			$functions->{$seedid} =~ s/\s*#.+//;
	  			$fnarray = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$functions->{$seedid}))];
	  			if ($combo eq join("\t",@{$fnarray})) {
	  				$UniqueRoleCombinations->{$combo}->{Sources}->[$index] .= "|seed";
	  			} else {
	  				if (defined($UniqueRoleCombinations->{$combo})) {
		  				$UniqueRoleCombinations->{$combo}->{Count}++;
		  				$index = @{$UniqueRoleCombinations->{$combo}->{SEEDGenes}};
		  				push(@{$UniqueRoleCombinations->{$combo}->{SEEDGenes}},$seedid);
		  				push(@{$UniqueRoleCombinations->{$combo}->{KBGenes}},$genes->[$j]);
		  				push(@{$UniqueRoleCombinations->{$combo}->{Sources}},"seed");
		  			} else {
		  				$index = 0;
		  				$UniqueRoleCombinations->{$combo} = {
		  					Count => 1,
		  					KBRep => $genes->[$j],
							SEEDRep => $seedid,
		  					KBGenes => [$genes->[$j]],
							SEEDGenes => [$seedid],
							Sources => ["seed"],
							Shared => {},
							Reactions => [],
							Roles => $fnarray,
							SS => []
		  				};
		  				$UniqueRoleCombinations->{$combo}->{RoleCount} = @{$fnarray};
		  				my $MetabolicCount = 0;
		  				my $SubsysCount = 0;
		  				my $UnannotatedCount = 0;
		  				my $Neighbors = 0;
		  				for (my $k=0; $k < @{$fnarray}; $k++) {
		  					$UniqueRoleCombinations->{$combo}->{SS}->[$k] = [];
  							$UniqueRoleCombinations->{$combo}->{Reactions}->[$k] = [];
		  					my $searchname = Bio::KBase::ObjectAPI::utilities::convertRoleToSearchRole($fnarray->[$k]);
		  					if (defined($RoleToSubsys->{$searchname})) {
		  						$UniqueRoleCombinations->{$combo}->{SS}->[$k] = [keys(%{$RoleToSubsys->{$searchname}})];
		  					}
		  					if (defined($RoleToRxn->{$searchname})) {
		  						my $rxns = [keys(%{$RoleToRxn->{$searchname}})];
		  						for (my $m=0; $m < @{$rxns}; $m++) {
		  							$UniqueRoleCombinations->{$combo}->{Reactions}->[$k]->[$m] = $rxns->[$m].":".$eqns->{$rxns->[$m]};
		  							for (my $n=($m+1); $n < @{$rxns}; $n++) {
		  								if (defined($rxnneighbors->{$rxns->[$m]}->{$rxns->[$n]})) {
		  									foreach my $cpd (keys(%{$rxnneighbors->{$rxns->[$m]}->{$rxns->[$n]}})) {
		 										$Neighbors++;
		 										my $cpdarr = [split(/_/,$cpd)];
		  										$UniqueRoleCombinations->{$combo}->{Shared}->{$cpdnames->{$cpdarr->[0]}} = 1;	
		  									}
		  								}
		  							}
		  						}
		  						$MetabolicCount++;
		  					} elsif (defined($RoleToSubsys->{$searchname})) {
		  						$SubsysCount++;
		  					} else {
		  						$UnannotatedCount++;
		  					}
		  				}
		  				if ($MetabolicCount > 0) {
		  					if ($MetabolicCount > 1) {
		  						if ($Neighbors > 0) {
		  							push(@{$NMM},$UniqueRoleCombinations->{$combo});
		  						} else {
		  							push(@{$NNMM},$UniqueRoleCombinations->{$combo});
		  						}
		  					} elsif ($SubsysCount > 0) {
		  						push(@{$SMSS},$UniqueRoleCombinations->{$combo});
		  					} else {
		  						push(@{$SMNS},$UniqueRoleCombinations->{$combo});
		  					}
		  				} elsif ($SubsysCount > 0) {
		  					if ($SubsysCount > 1) {
		  						push(@{$NMTS},$UniqueRoleCombinations->{$combo});
		  					} else {
		  						push(@{$NMSS},$UniqueRoleCombinations->{$combo});
		  					}
		  				} else {
		  					push(@{$NMNS},$UniqueRoleCombinations->{$combo});
		  				}
		  			}
	  			}
  			}
  		}
  		$genes = [];
  		#$i = @{$data};	
	}
}

open ( my $nmmfile, ">", $directory."/NeighboringMM.txt");
print $nmmfile "KBase Genes\tSEED Genes\tRole count\tCount\tShared Reactants\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$NMM}; $i++) {
	my $item = $NMM->[$i];
	print $nmmfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t".$item->{Count}."\t\"".join("\n",keys(%{$item->{Shared}}));
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $nmmfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $nmmfile "\n".$item->{Reactions}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $nmmfile "\n".$item->{SS}->[$j]->[$k];
		}
	}
	print $nmmfile "\"\n";
}
close($nmmfile);

open ( my $nnmmfile, ">", $directory."/NonneighboringMM.txt");
print $nnmmfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$NNMM}; $i++) {
	my $item = $NNMM->[$i];
	print $nnmmfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $nnmmfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $nnmmfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $nnmmfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $nnmmfile "\"\n";
}
close($nnmmfile);

open ( my $smssfile, ">", $directory."/SingleMetabolicSubsystem.txt");
print $smssfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$SMSS}; $i++) {
	my $item = $SMSS->[$i];
	print $smssfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $smssfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $smssfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $smssfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $smssfile "\"\n";
}
close($smssfile);

open ( my $smnsfile, ">", $directory."/SingleMetabolicNoSS.txt");
print $smnsfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$SMNS}; $i++) {
	my $item = $SMNS->[$i];
	print $smnsfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $smnsfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $smnsfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $smnsfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $smnsfile "\"\n";
}
close($smnsfile);

open ( my $nmtsfile, ">", $directory."/NonmetabolicTwoSubsystem.txt");
print $nmtsfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$NMTS}; $i++) {
	my $item = $NMTS->[$i];
	print $nmtsfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $nmtsfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $nmtsfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $nmtsfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $nmtsfile "\"\n";
}
close($nmtsfile);

open ( my $nmssfile, ">", $directory."/NonmetabolicSingleSubsystem.txt");
print $nmssfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$NMSS}; $i++) {
	my $item = $NMSS->[$i];
	print $nmssfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $nmssfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $nmssfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $nmssfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $nmssfile "\"\n";
}
close($nmssfile);

open ( my $nmnsfile, ">", $directory."/NonmetabolicNoSubsystem.txt");
print $nmnsfile "KBase Representative\tSEED Representative\tRole count\tCount\t\"Roles\nSubsystems\nReactions\"\n";
for (my $i=0; $i < @{$NMNS}; $i++) {
	my $item = $NMNS->[$i];
	print $nmnsfile $item->{KBRep}."\t".$item->{SEEDRep}."\t".$item->{RoleCount}."\t\"".$item->{Count};
	for (my $j=0; $j < @{$item->{Roles}}; $j++) {
		print $nmnsfile "\"\t\"".$item->{Roles}->[$j];
		for (my $k=0;$k < @{$item->{SS}->[$j]}; $k++) {
			print $nmnsfile "\n".$item->{SS}->[$j]->[$k];
		}
		for (my $k=0;$k < @{$item->{Reactions}->[$j]}; $k++) {
			print $nmnsfile "\n".$item->{Reactions}->[$j]->[$k];
		}
	}
	print $nmnsfile "\"\n";
}
close($nmnsfile);