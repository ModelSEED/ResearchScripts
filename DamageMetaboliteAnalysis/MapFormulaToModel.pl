use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

#Saving reaction count for database compounds
my $cpd_hash = Bio::KBase::utilities::compound_hash();
foreach my $cpd (keys(%{$cpd_hash})) {
	$cpd_hash->{$cpd}->{reaction_count} = 0;
}
my $rxn_hash = Bio::KBase::utilities::reaction_hash();
foreach my $rxn (keys(%{$rxn_hash})) {
	if (defined($rxn_hash->{$rxn}->{compound_ids})) {
		for (my $i=0; $i < @{$rxn_hash->{$rxn}->{compound_ids}}; $i++) {
			$cpd_hash->{$rxn_hash->{$rxn}->{compound_ids}->[$i]}->{reaction_count}++;
		}
	}
}

#Parsing arguments
my $workspace = $ARGV[0];
my $currmdl = $ARGV[1];
my $previous = $ARGV[2];
my $prefix = $ARGV[3];
my $SEEDName = $ARGV[4];
my $fresh_seed = $ARGV[5];
my $previous_peaks = {};
my $previous_cpds = {};
my $previous_others = {};
my $inchi_hash = {};
my $smiles_hash = {};
#Loading previous peaks
if (defined($previous)) {
	if (-e "/Users/chenry/Dropbox/workspace/PNNLSFA/".$previous.".peaks") {
		my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/".$previous.".peaks");
		for (my $i=1; $i < @{$lines}; $i++) {
			my $array = [split(/\t/,$lines->[$i])];
			my $subarray = [split(/;/,$array->[3])];
			for (my $j=0; $j < @{$subarray}; $j++) {
				$previous_peaks->{$array->[0]}->{$subarray->[$j]} = 1;
				if ($subarray->[$j] =~ m/cpd\d\d\d\d\d_/) {
					$previous_cpds->{$subarray->[$j]}->{$array->[0]} = 1;
				} else {
					$previous_others->{$subarray->[$j]}->{$array->[0]} = 1;
				}
			}
		}
	}
	
	my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/".$previous.".compounds");
	for (my $i=1; $i < @{$lines}; $i++) {
		my $array = [split(/\t/,$lines->[$i])];
		$inchi_hash->{$array->[0]} = {
			smiles => $array->[1],
			id => $array->[2],
			source => $array->[3]
		};
		$smiles_hash->{$array->[1]}->{$array->[2]} = $inchi_hash->{$array->[0]};
	}
}

#Processing model
my $eligible_count = 0;
my $ms_compound_count = 0;
my $model = $impl->util_get_object(Bio::KBase::utilities::buildref($currmdl,$workspace));
my $cpds = $model->modelcompounds();
my $compound_count = @{$cpds};
my $smiles_list = [];
my $id_list = [];
my $translation = {};
my $masslist = [];
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->id();
	push(@{$masslist},$cpds->[$i]->compound()->mass());
	if (length($cpds->[$i]->smiles()) > 0) {
		my $exists = 0;
		if (length($cpds->[$i]->inchikey()) > 0 && defined($inchi_hash->{$cpds->[$i]->inchikey()})) {
			$exists = 1;
			if ($cpds->[$i]->id() ne $inchi_hash->{$cpds->[$i]->inchikey()}->{id}) {
				$cpds->[$i]->id($inchi_hash->{$cpds->[$i]->inchikey()}->{id});
			}
			if ($inchi_hash->{$cpds->[$i]->inchikey()}->{source} eq "DB") {
				$inchi_hash->{$cpds->[$i]->inchikey()}->{source} = $currmdl;
			}
		} elsif ($cpds->[$i]->id() !~ m/cpd\d\d\d\d\d_/ && defined($smiles_hash->{$cpds->[$i]->smiles()})) {
			my $selected;
			foreach my $newcpd (keys(%{$smiles_hash->{$cpds->[$i]->smiles()}})) {
				if ($newcpd =~ m/cpd\d\d\d\d\d_/) {
					my $cpddata = $cpd_hash->{$newcpd};
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
			}
			if (defined($selected)) {
				$exists = 1;
				if  ($cpds->[$i]->id() ne $selected->{id}) {
					$cpds->[$i]->id($selected->{id});
					$cpds->[$i]->inchikey($selected->{inchikey});
					if ($smiles_hash->{$cpds->[$i]->smiles()}->{$selected->{id}}->{source} eq "DB") {
						$smiles_hash->{$cpds->[$i]->smiles()}->{$selected->{id}}->{source} = $currmdl;
					}
				}
			}
		}
		#Replacing generic IDs with specific IDs
		if ($cpds->[$i]->id() =~ m/pkc/) {
			my $newid = $cpds->[$i]->id();
			$newid =~ s/pkc/$prefix/;
			$cpds->[$i]->id($newid);
		}
		#Adding compound to compound database
		if ($exists == 0) {
			$inchi_hash->{$cpds->[$i]->inchikey()} = {
				id => $cpds->[$i]->id(),
				smiles => $cpds->[$i]->smiles(),
				source => $currmdl
			};
			$smiles_hash->{$cpds->[$i]->smiles()}->{$cpds->[$i]->id()} = 1;
		}
		$eligible_count++;
		push(@{$smiles_list},$cpds->[$i]->smiles());
		push(@{$id_list},$id);
	}
	#Saving modified IDs in translation hash
	if ($id ne $cpds->[$i]->id()) {
		$translation->{$id} = $cpds->[$i]->id();
	}
	if ($cpds->[$i]->id() =~ m/cpd\d\d\d\d\d_/) {
		$ms_compound_count++;
	}
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/".$currmdl.".mass",$masslist);

#Updating compound refs in reactions 
my $rxns = $model->modelreactions();
my $reaction_count = @{$rxns};
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/$prefix/;
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

#Saving modified model
my $wsmeta = $impl->util_save_object($model,$workspace."/".$currmdl.".fix");

#Loading or created seed model
my $seed_model;
if (defined($fresh_seed) && $fresh_seed == 1) {
	$seed_model = Bio::KBase::ObjectAPI::KBaseFBA::FBAModel->new({
		id => $SEEDName,
		source => "Pickaxe",
		source_id => $SEEDName,
		type => "CheminformaticsSeedModel",
		name => $SEEDName,
		template_ref => "NewKBaseModelTemplates/GramNegModelTemplateV2",
		template_refs => ["NewKBaseModelTemplates/GramNegModelTemplateV2"]
	});
} else {
	$seed_model = $impl->util_get_object(Bio::KBase::utilities::buildref($SEEDName,$workspace));
}

#Loading metabolomics and populating formula hash
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/MetabolomicsData.tsv");
my $formula_hash = {};
my $peak_hash = {};
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	$peak_hash->{$array->[0]} = [$array->[1],$array->[2]];
	if (length($array->[2]) > 0) {
		$formula_hash->{$array->[2]}->{$array->[0]} = $array->[1];
	}
}

#Converting model smiles into formula
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/temp/smiles.smi",$smiles_list);
system('obabel /Users/chenry/temp/smiles.smi -oinchi -T /nochg/formula > /Users/chenry/temp/inchi');
$lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/temp/inchi");

#Computing hits in the current model
my $peak_hits = {};
my $cpd_hits = {};
my $other_hits = {};
for (my $i=0; $i < @{$lines}; $i++) {
	if ($lines->[$i] =~ m/InChI.1S\/([^\/]+)\//) {
		my $formula = $1;
		if (defined($formula_hash->{$formula})) {
			foreach my $peakid (keys(%{$formula_hash->{$formula}})) {
				#Translating ID if a transaltion exists
				my $cpdid = $id_list->[$i];
				if (defined($translation->{$cpdid})) {
					$cpdid = $translation->{$cpdid};
				}
				#Adding hit to seed model
				if (!defined($seed_model->getObject("modelcompounds",$cpdid))) {
					$seed_model->add("modelcompounds",$model->getObject("modelcompounds",$cpdid)->cloneObject());
				}
				#Saving hits
				$peak_hits->{$peakid}->{$cpdid} = 1;
				if ($cpdid =~ m/cpd\d\d\d\d\d_/) {
					$cpd_hits->{$cpdid}->{$peakid} = 1;
				} else {
					$other_hits->{$cpdid}->{$peakid} = 1;
				}
			}
		}
	}
}
#$cpds = $seed_model->modelcompounds();
#for 

#Counting hits from the current model
my $model_peaks = keys(%{$peak_hits});
my $model_cpd = keys(%{$cpd_hits});
my $model_other = keys(%{$other_hits});
#Counting the peaks that are new in the current model
my $new_peaks = 0;
my $new_cpd = 0;
my $new_other = 0;
#Merging current model hits with previous hits
my $other_model_peaks = 0;
my $other_new_peaks = 0;
foreach my $peakid (keys(%{$peak_hits})) {
	my $new = 0;
	if (!defined($previous_peaks->{$peakid})) {
		$new = 1;
		$new_peaks++;
	}
	foreach my $cpdid (keys(%{$peak_hits->{$peakid}})) {
		$previous_peaks->{$peakid}->{$cpdid} = 1;
		if ($cpdid =~ m/cpd\d\d\d\d\d_/) {
			if (!defined($previous_cpds->{$cpdid})) {
				$new_cpd++;
			}
			$previous_cpds->{$cpdid}->{$peakid} = 1;
		} else {
			if (!defined($previous_others->{$cpdid})) {
				$new_other++;
			}
			$previous_others->{$cpdid}->{$peakid} = 1;
		}
	}
	my $ms_cpd = 0;
	foreach my $cpdid (keys(%{$previous_peaks->{$peakid}})) {
		if ($cpdid =~ m/cpd\d\d\d\d\d_/) {
			$ms_cpd = 1;
			last;
		}
	}
	if ($ms_cpd == 0) {
		$other_model_peaks++;
		if ($new == 1) {
			$other_new_peaks++;
		}
	}
}
my $total_peaks = keys(%{$previous_peaks});
my $total_cpd = keys(%{$previous_cpds});
my $total_other = keys(%{$previous_others});

#Reporting stats from current run
print "Total reactions:".$reaction_count."\n";
print "Total compounds:".$compound_count."\n";
print "Total MS compounds:".$ms_compound_count."\n";
print "Eligible compounds:".$eligible_count."\n";
print "Total peaks:".$total_peaks."\n";
print "Total cpd hits:".$total_cpd."\n";
print "Total other hits:".$total_other."\n";
print "Model peaks:".$model_peaks."\n";
print "No MS model peaks:".$other_model_peaks."\n";
print "Model cpd hits:".$model_cpd."\n";
print "Model other hits:".$model_other."\n";
print "New peaks:".$new_peaks."\n";
print "No MS new peaks:".$other_new_peaks."\n";
print "New cpd hits:".$new_cpd."\n";
print "New other hits:".$new_other."\n";

#Saving new hits file
my $output = ["id\tmass\tformula\tmsids"];
foreach my $peakid (keys(%{$previous_peaks})) {
	push(@{$output},$peakid."\t".$peak_hash->{$peakid}->[0]."\t".$peak_hash->{$peakid}->[1]."\t".join(";",keys(%{$previous_peaks->{$peakid}})));	
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/".$currmdl.".peaks",$output);

#Saving new compounds file
$output = ["inchikey\tsmiles\tid\tsource"];
foreach my $inchikey (keys(%{$inchi_hash})) {
	push(@{$output},$inchikey."\t".$inchi_hash->{$inchikey}->{smiles}."\t".$inchi_hash->{$inchikey}->{id}."\t".$inchi_hash->{$inchikey}->{source});
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/".$currmdl.".compounds",$output);

#Saving SEED model
if ($SEEDName ne "none") {
	$wsmeta = $impl->util_save_object($seed_model,$workspace."/".$SEEDName);
}