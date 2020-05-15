use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "46377";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("Cellvibrio_japonicus_mdl",$workspace));
my $dam_model = $impl->util_get_object(Bio::KBase::utilities::buildref("Spontaneous_Cellvibrio_expansion",$workspace));
my $rep_model = $impl->util_get_object(Bio::KBase::utilities::buildref("Enzyme_Cellvibrio_expansion",$workspace));

my $cpd_hash = {};
my $rejected_hash = {};
my $other_hash = {};
my $peak_string = "";
my $all_smiles = {};

#TODO: Compute average IDs per peak; Report stats from hashes

my $catagories = {
	db_hit => 0,
	model_hit => 0,
	spontaneous_hit => 0,
	enzyme_hit => 0,
	f_db_hit => 0,
	f_model_hit => 0,
	f_spontaneous_hit => 0,
	f_enzyme_hit => 0
};

my $model_cpd_ids = {};
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$model_cpd_ids->{$cpds->[$i]->id()} = 1;
	$all_smiles->{$cpds->[$i]->smiles()} = 1;
}
$cpds = $dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	#$all_smiles->{$cpds->[$i]->smiles()} = 1;
}
$cpds = $rep_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	#$all_smiles->{$cpds->[$i]->smiles()} = 1;
}
my $cpd_hash = Bio::KBase::utilities::compound_hash();
foreach my $cpd (keys(%{$cpd_hash})) {
	if (defined($cpd_hash->{$cpd}->{smiles})) {
		#$all_smiles->{$cpd_hash->{$cpd}->{smiles}} = 1;
	}
}
my $filearray = [];
my $count = 0;
foreach my $smiles (keys(%{$all_smiles})) {
	if (length($smiles) > 0) {
		push(@{$filearray},$smiles);
	}
}
#Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/smiles.smi",$filearray);
system('obabel /Users/chenry/Dropbox/workspace/PNNLSFA/smiles.smi -oinchi -T /nochg/formula > /Users/chenry/Dropbox/workspace/PNNLSFA/temp');
exit();

#Loading metabolomic peak inchi into hash
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/Cjaponicus_peaks.tsv");
my $other_lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/Cjaponicus_fulldata_peaks.tsv");
my $peak_hash;
my $peak_ids = [];
my $formula_hash = {};
for (my $i=1; $i < @{$other_lines}; $i++) {
	my $array = [split(/[\t]/,$other_lines->[$i])];
	push(@{$peak_ids},[$array->[0],"peak.".$i]);
	$peak_hash->{$array->[0]} = {
		id => "peak.".$i,
		isotopic => $array->[1],
		mf => $array->[2],
		error => $array->[3],
		ec => $array->[4],
		defect => $array->[5],
		kendrick_m => $array->[6],
		kendrick_defect => $array->[7],
		dbe => $array->[8],
		dbe_chnos => $array->[9],
		nso => $array->[10],
		c => $array->[11],
		o => $array->[12],
		h => $array->[13],
		n => $array->[14],
		s => $array->[15],
		p => $array->[16],
		aromaticity => $array->[17],
		cho => $array->[18],
		comment => $array->[19],
		TCP => $array->[20],
		TFM => $array->[21],
		TS => $array->[22],
		CCP => $array->[23],
		CFM => $array->[24],
		CS => $array->[25],
		HCP => $array->[26],
		HFM => $array->[27],
		HS => $array->[28]
	};
	if (defined($peak_hash->{$array->[0]}->{mf}) && length($peak_hash->{$array->[0]}->{mf}) > 0) {
		my $carbon = 0;
		if ($peak_hash->{$array->[0]}->{mf} =~ m/13C/) {
			$carbon = 1;
			$peak_hash->{$array->[0]}->{mf} =~ s/13C//;
		}
		if ($peak_hash->{$array->[0]}->{mf} =~ m/(.*)C(\d+)([A-Z]*.*)/) {
			$carbon += $2;
			$peak_hash->{$array->[0]}->{mf} = $1."C".$carbon.$3;
		} elsif ($peak_hash->{$array->[0]}->{mf} =~ m/(.*)C([A-Z]*.*)/) {
			$carbon += 1;
			$peak_hash->{$array->[0]}->{mf} = $1."C".$carbon.$3;
		} elsif ($carbon == 1) {
			$peak_hash->{$array->[0]}->{mf} = "C".$peak_hash->{$array->[0]}->{mf};
		}
		$formula_hash->{$peak_hash->{$array->[0]}->{mf}}->{$array->[0]} = 1;
	}
}
my $smiles_hash = {};
my $peak_hits_hash = {};
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/[\t]/,$lines->[$i])];
	$smiles_hash->{$array->[4]}->{$array->[2]} = 1;
	$peak_hash->{$array->[2]}->{smiles}->{$array->[4]} = 1;
}
#Building structure hashes
my $hash = {};
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $smiles = $cpds->[$i]->smiles();
	my $formula = &smiles_to_formula($smiles);
	if (defined($smiles_hash->{$smiles})) {
		foreach my $peak (keys(%{$smiles_hash->{$smiles}})) {
			$peak_hash->{$peak}->{ids}->{$cpds->[$i]->id()} = 1;
			$catagories->{model_hit}++;
			$cpd_hash->{peaks}->{$peak} = 1;
		}
		$cpd_hash->{ids}->{$cpds->[$i]->id()} = 1;
	}
	if (defined($formula_hash->{$formula})) {
		foreach my $peak (keys(%{$formula_hash->{$formula}})) {
			$peak_hash->{$peak}->{f_ids}->{$cpds->[$i]->id()} = 1;
			$catagories->{f_model_hit}++;
			$cpd_hash->{f_peaks}->{$peak} = 1;
		}
		$cpd_hash->{f_ids}->{$cpds->[$i]->id()} = 1;
	}
	$hash->{model}->{smiles}->{$smiles}->{$cpds->[$i]->id()} = $cpds->[$i];
	$hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()} = $cpds->[$i];
}

foreach my $cpd (keys(%{$cpd_hash})) {
	$cpd_hash->{$cpd}->{reaction_count} = 0;
	if (defined($cpd_hash->{$cpd}->{smiles})) {
		my $formula = &smiles_to_formula($cpd_hash->{$cpd}->{smiles});
		if (defined($smiles_hash->{$cpd_hash->{$cpd}->{smiles}})) {
			foreach my $peak (keys(%{$smiles_hash->{$cpd_hash->{$cpd}->{smiles}}})) {
				$peak_hash->{$peak}->{ids}->{$cpd_hash->{$cpd}->{id}."_c0"} = 1;
				$catagories->{db_hit}++;
				$cpd_hash->{peaks}->{$peak} = 1;
			}
			$cpd_hash->{ids}->{$cpd_hash->{$cpd}->{id}} = 1;
		}
		if (defined($formula_hash->{$formula})) {
			foreach my $peak (keys(%{$formula_hash->{$formula}})) {
				$peak_hash->{$peak}->{f_ids}->{$cpd_hash->{$cpd}->{id}."_c0"} = 1;
				$catagories->{f_db_hit}++;
				$cpd_hash->{f_peaks}->{$peak} = 1;
			}
			$cpd_hash->{f_ids}->{$cpd_hash->{$cpd}->{id}} = 1;
		}
		$hash->{db}->{smiles}->{$cpd_hash->{$cpd}->{smiles}}->{$cpd_hash->{$cpd}->{id}} = $cpd_hash->{$cpd};
		$hash->{db}->{inchikey}->{$cpd_hash->{$cpd}->{inchikey}} = $cpd_hash->{$cpd};
	}
}
my $rxn_hash = Bio::KBase::utilities::reaction_hash();
foreach my $rxn (keys(%{$rxn_hash})) {
	if (defined($rxn_hash->{$rxn}->{compound_ids})) {
		for (my $i=0; $i < @{$rxn_hash->{$rxn}->{compound_ids}}; $i++) {
			$cpd_hash->{$rxn_hash->{$rxn}->{compound_ids}->[$i]}->{reaction_count}++;
		}
	}
}

#Renaming IDs and combining the damage and repair model
my $translation = {};
$cpds = $dam_model->modelcompounds();
my $model_cpd_hash;
for (my $i=0; $i < @{$cpds}; $i++) {
	my $smiles = $cpds->[$i]->smiles();
	my $formula = &smiles_to_formula($smiles);
	if ($cpds->[$i]->id() =~ m/pkc/ && $cpds->[$i]->id() !~ m/_e0/ && $cpds->[$i]->inchikey() ne "none") {
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		if (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
			print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":model\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
		} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
			print $cpds->[$i]->id().":missed full match:".$cpddata->{id}.":".$cpddata->{name}.":db\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->{id}."_c0";
			$cpds->[$i]->id($cpddata->{id}."_c0");
		} elsif (defined($hash->{model}->{smiles}->{$smiles})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{model}->{smiles}->{$smiles}})) {
				my $cpddata = $hash->{model}->{smiles}->{$smiles}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
				print $cpds->[$i]->id().":missed base match:".$cpddata->id().":".$cpddata->name().":model\n";
			}
			print $cpds->[$i]->id().":selected missed base match:".$selected->id().":".$selected->name().":model\n";
			$model_cpd_hash->{$cpds->[$i]->id()} = $cpds->[$i];
			$translation->{$cpds->[$i]->id()} = $selected->id();
			$cpds->[$i]->id($selected->id());
		} elsif (defined($hash->{db}->{smiles}->{$smiles})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{db}->{smiles}->{$smiles}})) {
				my $cpddata = $hash->{db}->{smiles}->{$smiles}->{$newcpd};
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
				print $cpds->[$i]->id().":missed base match:".$cpddata->{id}.":".$cpddata->{name}.":".$cpddata->{reaction_count}.":db\n";
			}
			print $cpds->[$i]->id().":selected missed base match:".$selected->{id}.":".$selected->{name}.":".$selected->{reaction_count}.":db\n";
			$translation->{$cpds->[$i]->id()} = $selected->{id}."_c0";
			$cpds->[$i]->id($selected->{id}."_c0");
		}
	}
	if (defined($smiles_hash->{$smiles})) {
		foreach my $peak (keys(%{$smiles_hash->{$smiles}})) {
			$peak_hash->{$peak}->{ids}->{$cpds->[$i]->id()} = 1;
			$catagories->{spontaneous_hit}++;
			if ($cpds->[$i]->id() =~ m/cpd\d+/) {
				$cpd_hash->{peaks}->{$peak} = 1;
			} else {
				$other_hash->{peaks}->{$peak} = 1;
			}
		}
		if ($cpds->[$i]->id() =~ m/cpd\d+/) {
			$cpd_hash->{ids}->{$cpds->[$i]->id()} = 1;
		} else {
			$other_hash->{ids}->{$cpds->[$i]->id()} = 1;
		}
	}
	if (defined($formula_hash->{$formula})) {
		foreach my $peak (keys(%{$formula_hash->{$formula}})) {
			$peak_hash->{$peak}->{f_ids}->{$cpds->[$i]->id()} = 1;
			$catagories->{f_spontaneous_hit}++;
			if ($cpds->[$i]->id() =~ m/cpd\d+/) {
				$cpd_hash->{f_peaks}->{$peak} = 1;
			} else {
				$other_hash->{f_peaks}->{$peak} = 1;
			}
		}
		if ($cpds->[$i]->id() =~ m/cpd\d+/) {
			$cpd_hash->{f_ids}->{$cpds->[$i]->id()} = 1;
		} else {
			$other_hash->{f_ids}->{$cpds->[$i]->id()} = 1;
		}
	} else {#if ($cpds->[$i]->id() !~ m/cpd\d+/) {
		$rejected_hash->{$cpds->[$i]->id()} = 1;
	}
	$hash->{spont}->{inchikey}->{$cpds->[$i]->inchikey()} = $cpds->[$i];
	$hash->{spont}->{smiles}->{$smiles}->{$cpds->[$i]->id()} = $cpds->[$i];
	if ($cpds->[$i]->id() =~ m/pkc/) {
		my $id = $cpds->[$i]->id();
		$id =~ s/pkc/spont/;
		$translation->{$cpds->[$i]->id()} = $id;
		$cpds->[$i]->id($id);
	}
}

my $rxns = $dam_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/spont/;
		$rxns->[$i]->id($id);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		if ($rgts->[$j]->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$rgts->[$j]->modelcompound_ref($1.$translation->{$2});
			}
		}
	}
}

#Cloaning the damage model to become our new combined model
my $combined = $dam_model->cloneObject();
$combined->parent($dam_model->parent());
#Removing all transporters and extracellular metablites
$rxns = $combined->modelreactions();
foreach my $rxn (@{$rxns}) {
	if ($rxn->id() =~ m/transporter_c0/) {
		$combined->remove("modelreactions",$rxn);
	} else {
		my $rgts = $rxn->modelReactionReagents();
		my $remove = 0;
		for (my $j=0; $j < @{$rgts}; $j++) {
			if (defined($rejected_hash->{$rgts->[$j]->modelcompound()->id()})) {
				$remove = 1;
				last;
			}
		}
		if ($remove == 1) {
			$combined->remove("modelreactions",$rxn);
		}
	}
	$rxn->modelReactionProteins([]);
}
$cpds = $combined->modelcompounds();
foreach my $cpd (@{$cpds}) {
	if ($cpd->id() =~ m/_e0/ || defined($rejected_hash->{$cpd->id()})) {
		$combined->remove("modelcompounds",$cpd);
	}
}

$translation = {};
my $cpds = $rep_model->modelcompounds();
my $newcpdlist;
for (my $i=0; $i < @{$cpds}; $i++) {	
	my $smiles = $cpds->[$i]->smiles();
	my $formula = &smiles_to_formula($smiles);
	if ($cpds->[$i]->id() =~ m/pkc/ && $cpds->[$i]->id() !~ m/_e0/ && $cpds->[$i]->inchikey() ne "none") {
		if (defined($hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()};
			print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":damage\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
		} elsif (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
			print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":model\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
		} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
			print $cpds->[$i]->id().":missed full match:".$cpddata->{id}.":".$cpddata->{name}.":db\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->{id}."_c0";
			$cpds->[$i]->id($cpddata->{id}."_c0");
		} elsif (defined($hash->{model}->{smiles}->{$smiles})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{model}->{smiles}->{$smiles}})) {
				my $cpddata = $hash->{model}->{smiles}->{$smiles}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
				print $cpds->[$i]->id().":missed base match:".$cpddata->id().":".$cpddata->name().":model\n";
			}
			print $cpds->[$i]->id().":selected missed base match:".$selected->id().":".$selected->name().":model\n";
			$translation->{$cpds->[$i]->id()} = $selected->id();
			$cpds->[$i]->id($selected->id());
		} elsif (defined($hash->{db}->{smiles}->{$smiles})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{db}->{smiles}->{$smiles}})) {
				my $cpddata = $hash->{db}->{smiles}->{$smiles}->{$newcpd};
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
				print $cpds->[$i]->id().":missed base match:".$cpddata->{id}.":".$cpddata->{name}.":".$cpddata->{reaction_count}.":db\n";
			}
			print $cpds->[$i]->id().":selected missed base match:".$selected->{id}.":".$selected->{name}.":".$selected->{reaction_count}.":db\n";
			$translation->{$cpds->[$i]->id()} = $selected->{id}."_c0";
			$cpds->[$i]->id($selected->{id}."_c0");
		}
	}
	my $newcpd = $combined->getObject("modelcompounds",$cpds->[$i]->id());
	if (!defined($newcpd)) {	
		if ($cpds->[$i]->id() =~ m/pkc/) {
			my $id = $cpds->[$i]->id();
			$id =~ s/pkc/enz/;
			$translation->{$cpds->[$i]->id()} = $id;
			$cpds->[$i]->id($id);
		}
		if (defined($smiles_hash->{$smiles})) {
			if (defined($smiles_hash->{$smiles})) {
				foreach my $peak (keys(%{$smiles_hash->{$smiles}})) {
					$peak_hash->{$peak}->{ids}->{$cpds->[$i]->id()} = 1;
					$catagories->{enzyme_hit}++;
					if ($cpds->[$i]->id() =~ m/cpd\d+/) {
						$cpd_hash->{peaks}->{$peak} = 1;
					} else {
						$other_hash->{peaks}->{$peak} = 1;
					}
				}
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpd_hash->{ids}->{$cpds->[$i]->id()} = 1;
				} else {
					$other_hash->{ids}->{$cpds->[$i]->id()} = 1;
				}
			}
			#$combined->add("modelcompounds",$cpds->[$i]->cloneObject());
		}
		if (defined($formula_hash->{$formula})) {
			if (defined($formula_hash->{$formula})) {
				foreach my $peak (keys(%{$formula_hash->{$formula}})) {
					$peak_hash->{$peak}->{f_ids}->{$cpds->[$i]->id()} = 1;
					$catagories->{f_enzyme_hit}++;
					if ($cpds->[$i]->id() =~ m/cpd\d+/) {
						$cpd_hash->{f_peaks}->{$peak} = 1;
					} else {
						$other_hash->{f_peaks}->{$peak} = 1;
					}
				}
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpd_hash->{f_ids}->{$cpds->[$i]->id()} = 1;
				} else {
					$other_hash->{f_ids}->{$cpds->[$i]->id()} = 1;
				}
			}
			$combined->add("modelcompounds",$cpds->[$i]->cloneObject());
		} else {
			$rejected_hash->{$cpds->[$i]->id()} = 1;
		}
	}
}

$cpds = $combined->modelcompounds();
my $combined_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	$combined_hash->{$cpds->[$i]->id()} = 1;
}

for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();	
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/enz/;
		$rxns->[$i]->id($id);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		if ($rgts->[$j]->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$rgts->[$j]->modelcompound_ref($1.$translation->{$2});
			}
		}
	}
	if ($rxns->[$i]->id() !~ m/transporter_c0/) {
		my $rgts = $rxns->[$i]->modelReactionReagents();
		my $add = 1;
		for (my $j=0; $j < @{$rgts}; $j++) {
			if (!defined($combined_hash->{$rgts->[$j]->modelcompound()->id()})) {
				$add = 0;
				last;
			}
		}
		if ($add == 1) {
			$combined->add("modelreactions",$rxns->[$i]->cloneObject());
		}
	}
}

#my $wsmeta = $impl->util_save_object($dam_model,"29280/Spontaneous_Cellvibrio_expansion.fix");
#$wsmeta = $impl->util_save_object($rep_model,"29280/Enzyme_Cellvibrio_expansion.fix");

#Add drain reactions for every compound not in the original model
$cpds = $combined->modelcompounds();
my $combined_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	$combined_hash->{$cpds->[$i]->id()} = 1;
	if (!defined($model_cpd_ids->{$cpds->[$i]->id()})) {
		my $id = $cpds->[$i]->id();
		$id =~ s/_[a-z]\d$//;
		$combined->add("modelreactions",{
			id => $id."-drain_c0",
			name => $cpds->[$i]->name()." drain",
			direction => ">",
			reaction_ref => "~/modelreactions/id/rxn00000",
			modelcompartment_ref => "~/modelcompartments/id/c0",
			modelReactionReagents => [{
				modelcompound_ref => "~/modelcompounds/id/".$cpds->[$i]->id(),
				coefficient => -1
			}]
		});
	}
}
$rxns = $combined->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		my $ref = $rgts->[$j]->modelcompound_ref();
		$ref =~ s/_e0/_c0/;
		$rgts->[$j]->modelcompound_ref($ref);
	}
	$rxns->[$i]->modelReactionProteins([]);
}

#Save the combined model
my $wsmeta = $impl->util_save_object($combined,"46377/CombinedExpansionModel2");

my $cpdids = keys(%{$cpd_hash->{ids}});
my $otherids = keys(%{$other_hash->{ids}});
my $cpdpeaks = keys(%{$cpd_hash->{peaks}});
my $otherpeaks = keys(%{$other_hash->{peaks}});

my $totalpeaks = 0;
my $average_peaks = 0;
foreach my $peak (keys(%{$peak_hash})) {
	if (defined($peak_hash->{$peak}->{ids})) {
		if (length($peak_string) > 0) {
			$peak_string .= ";";
		}
		$peak_string .= $peak_hash->{$peak}->{id}.":1:".join(":",keys(%{$peak_hash->{$peak}->{ids}}));
		my $count = keys(%{$peak_hash->{$peak}->{ids}});
		$totalpeaks++;
		$average_peaks += $count;
	}
}
$average_peaks = $average_peaks/$totalpeaks;

print "Extracellular:".$catagories->{extracellular}."\n";
print "Database hit:".$catagories->{db_hit}."\n";
print "Enzyme metabolite hit:".$catagories->{enzyme_hit}."\n";
print "Spontaneous metabolite hit:".$catagories->{spontaneous_hit}."\n";
print "Model metabolite hit:".$catagories->{model_hit}."\n";
print "Compound peaks:".$cpdpeaks."\n";
print "Other peaks:".$otherpeaks."\n";
print "Compound IDs:".$cpdids."\n";
print "Other IDs:".$otherids."\n";
print "All IDs:".($cpdids+$otherids)."\n";
print "All peaks:".$totalpeaks."\n";
print "Average hits per peak:".$average_peaks."\n";

$cpdids = keys(%{$cpd_hash->{f_ids}});
$otherids = keys(%{$other_hash->{f_ids}});
$cpdpeaks = keys(%{$cpd_hash->{f_peaks}});
$otherpeaks = keys(%{$other_hash->{f_peaks}});

my $totalpeaks = 0;
my $average_peaks = 0;
foreach my $peak (keys(%{$peak_hash})) {
	if (defined($peak_hash->{$peak}->{f_ids})) {
		if (length($peak_string) > 0) {
			$peak_string .= ";";
		}
		$peak_string .= $peak_hash->{$peak}->{id}.":1:".join(":",keys(%{$peak_hash->{$peak}->{f_ids}}));
		my $count = keys(%{$peak_hash->{$peak}->{f_ids}});
		$totalpeaks++;
		$average_peaks += $count;
	}
}
$average_peaks = $average_peaks/$totalpeaks;

print "\nFORMULA DATA:\n";
print "Extracellular:".$catagories->{f_extracellular}."\n";
print "Database hit:".$catagories->{f_db_hit}."\n";
print "Enzyme metabolite hit:".$catagories->{f_enzyme_hit}."\n";
print "Spontaneous metabolite hit:".$catagories->{f_spontaneous_hit}."\n";
print "Model metabolite hit:".$catagories->{f_model_hit}."\n";
print "Compound peaks:".$cpdpeaks."\n";
print "Other peaks:".$otherpeaks."\n";
print "Compound IDs:".$cpdids."\n";
print "Other IDs:".$otherids."\n";
print "All IDs:".($cpdids+$otherids)."\n";
print "All peaks:".$totalpeaks."\n";
print "Average hits per peak:".$average_peaks."\n";

exit();
my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
	fbamodel_id => "Cellvibrio_japonicus_mdl",
	fbamodel_workspace => 46377,
	workspace => 46377,
	media_id => "Complete",
	media_workspace => "KBaseMedia",
	source_fbamodel_id => "CombinedExpansionModel",
	source_fbamodel_workspace => 46377,
	target_reaction => "bio1",
	fbamodel_output_id => "Cellvibrio_japonicus_mdl.fit",
	metabolite_peak_string => $peak_string
});

sub smiles_to_formula {
	my ($smiles) = @_;
	system('obabel -:"'.$smiles.'" -omolreport > /Users/chenry/Dropbox/workspace/PNNLSFA/temp');
	my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/temp");
	my $formula = "";
	for (my $i=0; $i < @{$lines}; $i++) {
		if ($lines->[$i] =~ m/FORMULA:\s(.+)/) {
			$formula = $1;
		}
	}
	return $formula;
}