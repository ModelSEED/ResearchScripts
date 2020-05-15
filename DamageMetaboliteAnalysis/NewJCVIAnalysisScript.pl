use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "29280";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn3",$workspace));
my $dam_model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn3_damage",$workspace));
my $rep_model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn3_repair",$workspace));
#Removing all transporters and extracellular metablites
my $cpds = $dam_model->modelcompounds();
foreach my $cpd (@{$cpds}) {
	if ($cpd->id() =~ m/_e0/) {
		$dam_model->remove("modelcompounds",$cpd);
	}
}
my $rxns = $dam_model->modelreactions();
foreach my $rxn (@{$rxns}) {
	if ($rxn->id() =~ m/transporter_c0/) {
		$dam_model->remove("modelreactions",$rxn);
	}
}
$cpds = $rep_model->modelcompounds();
foreach my $cpd (@{$cpds}) {
	if ($cpd->id() =~ m/_e0/) {
		$rep_model->remove("modelcompounds",$cpd);
	}
}
$rxns = $rep_model->modelreactions();
foreach my $rxn (@{$rxns}) {
	if ($rxn->id() =~ m/transporter_c0/) {
		$rep_model->remove("modelreactions",$rxn);
	}
}

#Loading metabolomic peak inchi into hash
#my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/Peaklist.txt");
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/FinalModeling/peak_list.txt");
my $peakstring = "";
my $peak_hash = {};
my $inchikey_hash = {};
my $metabolomics_hash = {};
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	$peak_hash->{$array->[0]} = {
		id => $array->[0],
		ave_rt => $array->[1],
		ave_mz => $array->[2],
		name => $array->[3],
		msi_level => $array->[4],
		sample_vs_media => $array->[5],
		adduct => $array->[6],
		formula => $array->[7],
		inchikey => $array->[8],
		pubchem => $array->[9],
		kegg => $array->[10],
		smiles => $array->[11],
		mass_error => $array->[12],
		platform => $array->[13]
	};
	my $keyarray = [split(/[_-]/,$peak_hash->{$array->[0]}->{inchikey})];
	$metabolomics_hash->{$keyarray->[0]}->{$array->[0]} = $peak_hash->{$array->[0]};
	$inchikey_hash->{$keyarray->[0]}->{$array->[0]} = $peak_hash->{$array->[0]};
}

#Building structure hashes
my $peakhits = {};
my $cpdhits = {};
my $otherhits = {};
my $totalcpds = {};
my $hash = {};
my $model_cpd_ids = {};
my $output_data = {};
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$model_cpd_ids->{$cpds->[$i]->id()} = 1;
	$totalcpds->{$cpds->[$i]->id()} = 1;
	if ($cpds->[$i]->id() =~ m/_c0$/) {
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				$metabolomics_hash->{$array->[0]}->{$peakid}->{model}->{$cpds->[$i]->id()} = 1;
				$peakhits->{all}->{$peakid}->{$cpds->[$i]->id()} = 1;
				$peakhits->{model}->{$peakid}->{$cpds->[$i]->id()} = 1;
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpdhits->{model}->{$cpds->[$i]->id()}->{$peakid} = 1;
					$cpdhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				} else {
					$otherhits->{model}->{$cpds->[$i]->id()}->{$peakid} = 1;
					$otherhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				}
			}
		}
		$hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()} = $cpds->[$i];
		$hash->{model}->{base}->{$array->[0]}->{$cpds->[$i]->id()} = $cpds->[$i];
	}
}
$output_data->{model_compounds} = @{$cpds};
$output_data->{model_reactions} = @{$model->modelreactions()};
$output_data->{total_model_peaks} = keys(%{$peakhits->{model}});
$output_data->{total_model_peak_compounds} = keys(%{$cpdhits->{model}})+keys(%{$otherhits->{model}});

my $cpd_hash = Bio::KBase::utilities::compound_hash();
foreach my $cpd (keys(%{$cpd_hash})) {
	$cpd_hash->{$cpd}->{reaction_count} = 0;
	if (defined($cpd_hash->{$cpd}->{inchikey})) {
		my $array = [split(/[_-]/,$cpd_hash->{$cpd}->{inchikey})];
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				$metabolomics_hash->{$array->[0]}->{$peakid}->{db}->{$cpd."_c0"} = 1;
				$peakhits->{db}->{$peakid}->{$cpd} = 1;
				$cpdhits->{db}->{$cpd."_c0"}->{$peakid} = 1;
			}
		}
		$hash->{db}->{inchikey}->{$cpd_hash->{$cpd}->{inchikey}} = $cpd_hash->{$cpd};
		$hash->{db}->{base}->{$array->[0]}->{$cpd} = $cpd_hash->{$cpd};
	}
}
$output_data->{total_db_peaks} = keys(%{$peakhits->{db}});
$output_data->{total_db_peak_compounds} = keys(%{$cpdhits->{db}});

my $rxn_hash = Bio::KBase::utilities::reaction_hash();
foreach my $rxn (keys(%{$rxn_hash})) {
	if (defined($rxn_hash->{$rxn}->{compound_ids})) {
		for (my $i=0; $i < @{$rxn_hash->{$rxn}->{compound_ids}}; $i++) {
			$cpd_hash->{$rxn_hash->{$rxn}->{compound_ids}->[$i]}->{reaction_count}++;
		}
	}
}

$cpds = $dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
	$hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()} = $cpds->[$i];
	$hash->{dam}->{base}->{$array->[0]}->{$cpds->[$i]->id()} = $cpds->[$i];
}

my $translation = {};
my $cpds = $dam_model->modelcompounds();
my $model_cpd_hash;
my $cpd_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	my $original_id = $cpds->[$i]->id();
	my $source = "damage";
	if (length($cpds->[$i]->inchikey()) > 0 && $cpds->[$i]->inchikey() ne "none") {
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		if (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
			$translation->{$original_id} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
			$source = "model";
		} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
			$translation->{$original_id} = $cpddata->{id}."_c0";
			$cpds->[$i]->id($cpddata->{id}."_c0");
			$source = "db";
		} elsif (defined($hash->{model}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{model}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{model}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
			}
			$translation->{$original_id} = $selected->id();
			$cpds->[$i]->id($selected->id());
			$source = "model";
		} elsif (defined($hash->{db}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{db}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{db}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $cpddata->{reaction_count} >= $selected->{reaction_count}) {
					if (!defined($selected)) {
						$selected = $cpddata;				
					} elsif ($cpddata->{reaction_count} == $selected->{reaction_count}) {
						if ($cpddata->{id} < $selected->{id}) {
							$selected = $cpddata;
						}
					} else {
						$selected = $cpddata;
					}
				}
			}
			$translation->{$original_id} = $selected->{id}."_c0";
			$cpds->[$i]->id($selected->{id}."_c0");
			$source = "db";
		}
		if ($cpds->[$i]->id() =~ m/pkc/) {
			my $id = $cpds->[$i]->id();
			$id =~ s/pkc/dmg/;
			$translation->{$cpds->[$i]->id()} = $id;
			$cpds->[$i]->id($id);
		}
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				if ($source ne "model") {
					$metabolomics_hash->{$array->[0]}->{$peakid}->{damage}->{$cpds->[$i]->id()} = 1;
					$peakhits->{damage}->{$peakid}->{$cpds->[$i]->id()} = 1;
					if ($cpds->[$i]->id() =~ m/cpd\d+/) {
						$cpdhits->{damage}->{$cpds->[$i]->id()}->{$peakid} = 1;
					} else {
						$otherhits->{damage}->{$cpds->[$i]->id()}->{$peakid} = 1;
					}
				}
				$peakhits->{all}->{$peakid}->{$cpds->[$i]->id()} = 1;
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpdhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				} else {
					if ($cpds->[$i]->name() =~ /pkc/) {
						$cpds->[$i]->name($peak_hash->{$peakid}->{name});
					}
					$otherhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				}
			}	
		}
	}
	if ($cpds->[$i]->id() =~ m/pkc/) {
		my $id = $cpds->[$i]->id();
		$id =~ s/pkc/dmg/;
		$translation->{$cpds->[$i]->id()} = $id;
		$cpds->[$i]->id($id);
	}
	if (!defined($model_cpd_hash->{$cpds->[$i]->id()})) {
		$model_cpd_hash->{$cpds->[$i]->id()} = $cpds->[$i];
	}
	$totalcpds->{$cpds->[$i]->id()} = 1;
	if ($cpds->[$i]->id() =~ m/(cpd\d+)/ && $cpds->[$i]->name() =~ m/^pkc/) {
		my $cpdid = $1;
		if (defined($cpd_hash->{$cpdid}->{name})) {
			$cpds->[$i]->name($cpd_hash->{$cpdid}->{name});
		}
	}
	if (defined($cpd_hash->{$cpds->[$i]->id()})) {
		$dam_model->remove("modelcompounds",$cpds->[$i]);
	}
	$cpd_hash->{$cpds->[$i]->id()} = 1;
}
$output_data->{total_damage_peaks} = keys(%{$peakhits->{damage}});
$output_data->{total_damage_peak_compounds} = keys(%{$cpdhits->{damage}});
$rxns = $dam_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/dmg/;
		$rxns->[$i]->id($id);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		my $modelcompound_ref = $rgts->[$j]->modelcompound_ref();
		$modelcompound_ref =~ s/_e0/_c0/;
		$rgts->[$j]->modelcompound_ref($modelcompound_ref);
		if ($rgts->[$j]->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$rgts->[$j]->modelcompound_ref($1.$translation->{$2});
			}
		}
	}
}
$output_data->{total_damage_rxn} = @{$rxns};
$output_data->{total_damage_cpd} = @{$cpds};
$output_data->{total_damage_and_model} = keys(%{$totalcpds});
$output_data->{total_damage_and_model_peaks} = keys(%{$peakhits->{all}});
$output_data->{total_damage_and_model_peak_compounds} = keys(%{$cpdhits->{all}}) + keys(%{$otherhits->{all}});

#Cloaning the damage model to become our new combined model
my $combined = $dam_model->cloneObject();
$combined->parent($dam_model->parent());

#Saving model
my $wsmeta = $impl->util_save_object($combined,"29280/CorrectedDamage");

#Processing repair model
$translation = {};
my $cpds = $rep_model->modelcompounds();
my $newcpdlist;
my $repcount = 0;
$cpd_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {	
	my $original_id = $cpds->[$i]->id();
	my $source = "repair";
	if (length($cpds->[$i]->inchikey()) > 0 && $cpds->[$i]->inchikey() ne "none") {
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		if (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
			$source = "model";
		} elsif (defined($hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()};
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
			$source = "damage";
		} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
			$translation->{$cpds->[$i]->id()} = $cpddata->{id}."_c0";
			$cpds->[$i]->id($cpddata->{id}."_c0");
		} elsif (defined($hash->{model}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{model}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{model}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
			}
			$translation->{$cpds->[$i]->id()} = $selected->id();
			$cpds->[$i]->id($selected->id());
			$source = "model";
		} elsif (defined($hash->{dam}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{dam}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{dam}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
			}
			$translation->{$cpds->[$i]->id()} = $selected->id();
			$cpds->[$i]->id($selected->id());
			$source = "damage";
		} elsif (defined($hash->{db}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{db}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{db}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $cpddata->{reaction_count} >= $selected->{reaction_count}) {
					if (!defined($selected)) {
						$selected = $cpddata;				
					} elsif ($cpddata->{reaction_count} == $selected->{reaction_count}) {
						if ($cpddata->{id} < $selected->{id}) {
							$selected = $cpddata;
						}
					} else {
						$selected = $cpddata;
					}
				}
			}
			$translation->{$cpds->[$i]->id()} = $selected->{id}."_c0";
			$cpds->[$i]->id($selected->{id}."_c0");
		}
		if ($cpds->[$i]->id() =~ m/pkc/) {
			my $id = $cpds->[$i]->id();
			$id =~ s/pkc/rep/;
			$translation->{$cpds->[$i]->id()} = $id;
			$cpds->[$i]->id($id);
		}
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				if ($source ne "model" && $source ne "damage") {
					$metabolomics_hash->{$array->[0]}->{$peakid}->{repair}->{$cpds->[$i]->id()} = 1;
					$peakhits->{repair}->{$peakid}->{$cpds->[$i]->id()} = 1;
					if ($cpds->[$i]->id() =~ m/cpd\d+/) {
						$cpdhits->{repair}->{$cpds->[$i]->id()}->{$peakid} = 1;
					} else {
						$otherhits->{repair}->{$cpds->[$i]->id()}->{$peakid} = 1;
					}
				}
				$peakhits->{all}->{$peakid}->{$cpds->[$i]->id()} = 1;
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpdhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				} else {
					if ($cpds->[$i]->name() =~ /pkc/) {
						$cpds->[$i]->name($peak_hash->{$peakid}->{name});
					}
					$otherhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				}
			}
		}
		my $newcpd = $combined->getObject("modelcompounds",$cpds->[$i]->id());
		if (!defined($newcpd) && ($cpds->[$i]->id() =~ m/cpd\d+/ || defined($metabolomics_hash->{$array->[0]}))) {
			$combined->add("modelcompounds",$cpds->[$i]->cloneObject());
		}
		$repcount++;
		$totalcpds->{$cpds->[$i]->id()} = 1;	
	}
	if ($cpds->[$i]->id() =~ m/(cpd\d+)/ && $cpds->[$i]->name() =~ m/^pkc/) {
		my $cpdid = $1;
		if (defined($cpd_hash->{$cpdid}->{name})) {
			$cpds->[$i]->name($cpd_hash->{$cpdid}->{name});
		}
	}
	if (defined($cpd_hash->{$cpds->[$i]->id()})) {
		$rep_model->remove("modelcompounds",$cpds->[$i]);
	}
	$cpd_hash->{$cpds->[$i]->id()} = 1;
}
$cpds = $combined->modelcompounds();
my $combined_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	$combined_hash->{$cpds->[$i]->id()} = 1;
}
my $rxns = $rep_model->modelreactions();
$output_data->{total_repair_rxn} = @{$rxns};
$output_data->{retained_repair_rxn} = 0;
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();	
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/rep/;
		$rxns->[$i]->id($id);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		my $modelcompound_ref = $rgts->[$j]->modelcompound_ref();
		$modelcompound_ref =~ s/_e0/_c0/;
		$rgts->[$j]->modelcompound_ref($modelcompound_ref);
		if ($rgts->[$j]->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$rgts->[$j]->modelcompound_ref($1.$translation->{$2});
			}
		}
	}
	$repcount++;
	my $rgts = $rxns->[$i]->modelReactionReagents();
	my $add = 1;
	for (my $j=0; $j < @{$rgts}; $j++) {
		if (!defined($combined_hash->{$rgts->[$j]->modelcompound()->id()})) {
			$add = 0;
			last;
		}
	}
	$rxns->[$i]->modelReactionProteins([]);
	if ($add == 1) {
		$output_data->{retained_repair_rxn}++;
		$combined->add("modelreactions",$rxns->[$i]->cloneObject());
	}
}
$output_data->{total_repair_cpd} = @{$cpds};
$output_data->{total_repair_and_damage_and_model} = keys(%{$totalcpds});
$output_data->{total_repair_damage_and_model_peaks} = keys(%{$peakhits->{all}});
$output_data->{total_repair_damage_and_model_peak_compounds} = keys(%{$cpdhits->{all}}) + keys(%{$otherhits->{all}});

#Saving model
$wsmeta = $impl->util_save_object($rep_model,"29280/CorrectedRepair");
#Saving model
$wsmeta = $impl->util_save_object($combined,"29280/Combined");

foreach my $label (keys(%{$output_data})) {
	print $label."\t".$output_data->{$label}."\n";
}

#foreach my $cpdid (keys(%{$cpdhits->{db}})) {
#	foreach my $peakid (keys(%{$cpdhits->{db}->{$cpdid}})) {
#		$cpdhits->{all}->{$cpdid}->{$peakid} = 1;
#		$peakhits->{all}->{$peakid}->{$cpdid} = 1;
#		$peakhits->{db}->{$peakid}->{$cpdid} = 1;
#	}
#}

##Add drain reactions for every compound not in the original model
#$cpds = $combined->modelcompounds();
#my $combined_hash = {};
#for (my $i=0; $i < @{$cpds}; $i++) {
#	$combined_hash->{$cpds->[$i]->id()} = 1;
#	if (!defined($model_cpd_ids->{$cpds->[$i]->id()})) {
#		my $id = $cpds->[$i]->id();
#		$id =~ s/_[a-z]\d$//;
#		$combined->add("modelreactions",{
#			id => $id."-drain_c0",
#			name => $cpds->[$i]->name()." drain",
#			direction => ">",
#			reaction_ref => "~/modelreactions/id/rxn00000",
#			modelcompartment_ref => "~/modelcompartments/id/c0",
#			modelReactionReagents => [{
#				modelcompound_ref => "~/modelcompounds/id/".$cpds->[$i]->id(),
#				coefficient => -1
#			}]
#		});
#	}
#}
#$rxns = $combined->modelreactions();
#for (my $i=0; $i < @{$rxns}; $i++) {
#	my $rgts = $rxns->[$i]->modelReactionReagents();
#	for (my $j=0; $j < @{$rgts}; $j++) {
#		my $ref = $rgts->[$j]->modelcompound_ref();
#		$ref =~ s/_e0/_c0/;
#		$rgts->[$j]->modelcompound_ref($ref);
#	}
#	$rxns->[$i]->modelReactionProteins([]);
#}

#Now adding original model to combined model'
my $fullcombined = $model->cloneObject();
$cpds = $combined->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (!defined($fullcombined->getObject("modelcompounds",$cpds->[$i]->id()))) {
		$fullcombined->add("modelcompounds",$cpds->[$i]->cloneObject());
	}
}
$rxns = $combined->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$fullcombined->add("modelreactions",$rxns->[$i]->cloneObject());
}

#Saving fully combined model
my $wsmeta = $impl->util_save_object($fullcombined,"29280/FullCombined.fix");

#Printing peak table
my $peak_tbl = ["Peak ID\tModel\tModelSEED\tDamage\tRepair"];

print "Model and damage only:\n";
foreach my $peakid (keys(%{$peakhits->{all}})) {
	print $peakid.":1";
	foreach my $cpdid (keys(%{$peakhits->{all}->{$peakid}})) {
		print ":".$cpdid;
	}
	print ";";
}
print "\n\nModel and damage and database:\n";
foreach my $peakid (keys(%{$peakhits->{db}})) {
	foreach my $cpdid (keys(%{$peakhits->{db}->{$peakid}})) {
		$peakhits->{all}->{$peakid}->{$cpdid} = 1;
	}
}
foreach my $peakid (keys(%{$peakhits->{all}})) {
	print $peakid.":1";
	foreach my $cpdid (keys(%{$peakhits->{all}->{$peakid}})) {
		print ":".$cpdid;
	}
	print ";";
}