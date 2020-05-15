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
my $old_dam_model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.damage-deleted-1572588332028",$workspace));
my $old_rep_model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.repair-deleted-1572588335426",$workspace));

my $peakhits = {};
my $cpdhits = {};
my $otherhits = {};
my $totalcpds = {};

my $model_cpd_ids = {};
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$model_cpd_ids->{$cpds->[$i]->id()} = 1;
}

#Loading metabolomic peak inchi into hash
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/Peaklist.txt");
my $peakstring = "";
my $peakids;
my $metabolomics_hash = {};
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/[\t\s]/,$lines->[$i])];
	my $subarray = [split(/[_-]/,$array->[0])];
	$metabolomics_hash->{$subarray->[0]}->{"peak.".($i+1)} = 1;
}

#Building structure hashes
my $hash = {};
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$totalcpds->{$cpds->[$i]->id()} = 1;
	my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
	if (defined($metabolomics_hash->{$array->[0]})) {
		foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
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
my $rxns = $model->modelreactions();

print "Total model cpd:".@{$cpds}."\n";
print "Total model rxn:".@{$rxns}."\n";
print "Model peaks:".keys(%{$peakhits->{model}})."\n";
print "Model seed compounds:".keys(%{$cpdhits->{model}})."\n";
print "Model other compounds:".keys(%{$otherhits->{model}})."\n";

my $cpd_hash = Bio::KBase::utilities::compound_hash();
foreach my $cpd (keys(%{$cpd_hash})) {
	$cpd_hash->{reaction_count} = 0;
	if (defined($cpd_hash->{$cpd}->{inchikey})) {
		my $array = [split(/[_-]/,$cpd_hash->{$cpd}->{inchikey})];
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				$cpdhits->{db}->{$cpd_hash->{$cpd}->{id}."_c0"}->{$peakid} = 1;
			}
		}
		$hash->{db}->{inchikey}->{$cpd_hash->{$cpd}->{inchikey}} = $cpd_hash->{$cpd};
		$hash->{db}->{base}->{$array->[0]}->{$cpd_hash->{$cpd}->{id}} = $cpd_hash->{$cpd};
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
$cpds = $dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
	$hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()} = $cpds->[$i];
	$hash->{dam}->{base}->{$array->[0]}->{$cpds->[$i]->id()} = $cpds->[$i];
}

#Saving the old damage and repair model IDs
my $original_cpd_hash = {};
my $original_cpd_id_hash = {};
my $cpds = $old_dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$original_cpd_id_hash->{damage}->{$cpds->[$i]->id()} = 0;
	$original_cpd_hash->{$cpds->[$i]->inchikey()}->{damage} = $cpds->[$i]->id();
}
$cpds = $old_rep_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	$original_cpd_id_hash->{repair}->{$cpds->[$i]->id()} = 0;
	$original_cpd_hash->{$cpds->[$i]->inchikey()}->{repair} = $cpds->[$i]->id();
}
my $original_rxn_hash = {};
my $original_rxn_id_hash = {};
my $rxns = $old_dam_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$original_rxn_id_hash->{damage}->{$rxns->[$i]->id()} = 0;
	$original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage} = $rxns->[$i]->id();
}
$rxns = $old_rep_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$original_rxn_id_hash->{repair}->{$rxns->[$i]->id()} = 0;
	$original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair} = $rxns->[$i]->id();
}

#Renaming IDs and combining the damage and repair model
my $rxns = $dam_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage})) {
		$rxns->[$i]->dblinks()->{olddmg} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage};
		$original_rxn_id_hash->{damage}->{$rxns->[$i]->id()} = 1;
	}
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair})) {
		$rxns->[$i]->dblinks()->{oldrep} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair};
		$original_rxn_id_hash->{repair}->{$rxns->[$i]->id()} = 1;
	}
}
my $translation = {};
my $cpds = $dam_model->modelcompounds();
my $model_cpd_hash;
for (my $i=0; $i < @{$cpds}; $i++) {
	if (length($cpds->[$i]->inchikey()) > 0 && $cpds->[$i]->inchikey() ne "none") {
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		if ($cpds->[$i]->id() =~ m/pkc/ && $cpds->[$i]->id() !~ m/_e0/ && $cpds->[$i]->inchikey() ne "none") {
			if (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
				my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
				#print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":model\n";
				$translation->{$cpds->[$i]->id()} = $cpddata->id();
				$cpds->[$i]->id($cpddata->id());
			} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
				my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
				#print $cpds->[$i]->id().":missed full match:".$cpddata->{id}.":".$cpddata->{name}.":db\n";
				$translation->{$cpds->[$i]->id()} = $cpddata->{id}."_c0";
				$cpds->[$i]->id($cpddata->{id}."_c0");
			} elsif (defined($hash->{model}->{base}->{$array->[0]})) {
				my $selected;
				foreach my $newcpd (keys(%{$hash->{model}->{base}->{$array->[0]}})) {
					my $cpddata = $hash->{model}->{base}->{$array->[0]}->{$newcpd};
					if (!defined($selected) || $selected->id() =~ m/_e0$/) {
						$selected = $cpddata;
					}
					#print $cpds->[$i]->id().":missed base match:".$cpddata->id().":".$cpddata->name().":model\n";
				}
				#print $cpds->[$i]->id().":selected missed base match:".$selected->id().":".$selected->name().":model\n";
				$translation->{$cpds->[$i]->id()} = $selected->id();
				$cpds->[$i]->id($selected->id());
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
					#print $cpds->[$i]->id().":missed base match:".$cpddata->{id}.":".$cpddata->{name}.":".$cpddata->{reaction_count}.":db\n";
				}
				#print $cpds->[$i]->id().":selected missed base match:".$selected->{id}.":".$selected->{name}.":".$selected->{reaction_count}.":db\n";
				$translation->{$cpds->[$i]->id()} = $selected->{id}."_c0";
				$cpds->[$i]->id($selected->{id}."_c0");
			}
		}
		if (defined($metabolomics_hash->{$array->[0]})) {
			foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
				$peakhits->{all}->{$peakid}->{$cpds->[$i]->id()} = 1;
				$peakhits->{dam}->{$peakid}->{$cpds->[$i]->id()} = 1;
				if ($cpds->[$i]->id() =~ m/cpd\d+/) {
					$cpdhits->{dam}->{$cpds->[$i]->id()}->{$peakid} = 1;
					$cpdhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
				} else {
					$otherhits->{dam}->{$cpds->[$i]->id()}->{$peakid} = 1;
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
	$model_cpd_hash->{$cpds->[$i]->id()} = $cpds->[$i];
	$totalcpds->{$cpds->[$i]->id()} = 1;
	if (defined($original_cpd_hash->{$cpds->[$i]->inchikey()}->{damage})) {
		$cpds->[$i]->dblinks()->{olddmg} = $original_cpd_hash->{$cpds->[$i]->inchikey()}->{damage};
		$original_cpd_id_hash->{damage}->{$rxns->[$i]->id()} = 1;
	}
	if (defined($original_cpd_hash->{$cpds->[$i]->inchikey()}->{repair})) {
		$cpds->[$i]->dblinks()->{oldrep} = $original_cpd_hash->{$cpds->[$i]->inchikey()}->{repair};
		$original_cpd_id_hash->{repair}->{$rxns->[$i]->id()} = 1;
	}
}
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/dmg/;
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
print "Total damage cpd:".keys(%{$model_cpd_hash})."\n";
print "Total damage rxn:".@{$rxns}."\n";
print "Total cpd:".keys(%{$totalcpds})."\n";
print "Damage peaks:".keys(%{$peakhits->{dam}})."\n";
print "Damage seed compounds:".keys(%{$cpdhits->{dam}})."\n";
print "Damage other compounds:".keys(%{$otherhits->{dam}})."\n";
print "Total peaks:".keys(%{$peakhits->{all}})."\n";
print "Total seed compounds:".keys(%{$cpdhits->{all}})."\n";
print "Total other compounds:".keys(%{$otherhits->{all}})."\n";

#Cloaning the damage model to become our new combined model
my $combined = $dam_model->cloneObject();
$combined->parent($dam_model->parent());
#Removing all transporters and extracellular metablites
$cpds = $combined->modelcompounds();
foreach my $cpd (@{$cpds}) {
	if ($cpd->id() =~ m/_e0/) {
		$combined->remove("modelcompounds",$cpd);
	}
}
$rxns = $combined->modelreactions();
foreach my $rxn (@{$rxns}) {
	if ($rxn->id() =~ m/transporter_c0/) {
		$combined->remove("modelreactions",$rxn);
	}
}

my $rxns = $rep_model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage})) {
		$rxns->[$i]->dblinks()->{olddmg} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage};
		$original_rxn_id_hash->{damage}->{$rxns->[$i]->id()} = 1;
	}
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair})) {
		$rxns->[$i]->dblinks()->{oldrep} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair};
		$original_rxn_id_hash->{repair}->{$rxns->[$i]->id()} = 1;
	}
}
$translation = {};
my $cpds = $rep_model->modelcompounds();
my $newcpdlist;
my $repcount = 0;
for (my $i=0; $i < @{$cpds}; $i++) {	
	my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
	if ($cpds->[$i]->id() =~ m/pkc/ && $cpds->[$i]->id() !~ m/_e0/ && $cpds->[$i]->inchikey() ne "none") {
		if (defined($hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{dam}->{inchikey}->{$cpds->[$i]->inchikey()};
			#print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":damage\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
		} elsif (defined($hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{model}->{inchikey}->{$cpds->[$i]->inchikey()};
			#print $cpds->[$i]->id().":missed full match:".$cpddata->id().":".$cpddata->name().":model\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->id();
			$cpds->[$i]->id($cpddata->id());
		} elsif (defined($hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()})) {
			my $cpddata = $hash->{db}->{inchikey}->{$cpds->[$i]->inchikey()};
			#print $cpds->[$i]->id().":missed full match:".$cpddata->{id}.":".$cpddata->{name}.":db\n";
			$translation->{$cpds->[$i]->id()} = $cpddata->{id}."_c0";
			$cpds->[$i]->id($cpddata->{id}."_c0");
		} elsif (defined($hash->{model}->{base}->{$array->[0]})) {
			my $selected;
			foreach my $newcpd (keys(%{$hash->{model}->{base}->{$array->[0]}})) {
				my $cpddata = $hash->{model}->{base}->{$array->[0]}->{$newcpd};
				if (!defined($selected) || $selected->id() =~ m/_e0$/) {
					$selected = $cpddata;
				}
				#print $cpds->[$i]->id().":missed base match:".$cpddata->id().":".$cpddata->name().":model\n";
			}
			#print $cpds->[$i]->id().":selected missed base match:".$selected->id().":".$selected->name().":model\n";
			$translation->{$cpds->[$i]->id()} = $selected->id();
			$cpds->[$i]->id($selected->id());
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
				#print $cpds->[$i]->id().":missed base match:".$cpddata->{id}.":".$cpddata->{name}.":".$cpddata->{reaction_count}.":db\n";
			}
			#print $cpds->[$i]->id().":selected missed base match:".$selected->{id}.":".$selected->{name}.":".$selected->{reaction_count}.":db\n";
			$translation->{$cpds->[$i]->id()} = $selected->{id}."_c0";
			$cpds->[$i]->id($selected->{id}."_c0");
		}
	}
	if ($cpds->[$i]->id() !~ m/_e0/) {
		$repcount++;
		$totalcpds->{$cpds->[$i]->id()} = 1;
	}
	my $newcpd = $combined->getObject("modelcompounds",$cpds->[$i]->id());
	if (!defined($newcpd)) {	
		if ($cpds->[$i]->id() =~ m/pkc/) {
			my $id = $cpds->[$i]->id();
			$id =~ s/pkc/rep/;
			$translation->{$cpds->[$i]->id()} = $id;
			$cpds->[$i]->id($id);
		}
		if ($cpds->[$i]->id() !~ m/_e0/ && ($cpds->[$i]->id() =~ m/cpd\d+/ || defined($metabolomics_hash->{$array->[0]}))) {
			if (defined($metabolomics_hash->{$array->[0]})) {
				foreach my $peakid (keys(%{$metabolomics_hash->{$array->[0]}})) {
					$peakhits->{all}->{$peakid}->{$cpds->[$i]->id()} = 1;
					$peakhits->{rep}->{$peakid}->{$cpds->[$i]->id()} = 1;
					if ($cpds->[$i]->id() =~ m/cpd\d+/) {
						$otherhits->{rep}->{$cpds->[$i]->id()}->{$peakid} = 1;
						$otherhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
					} elsif (defined($metabolomics_hash->{$array->[0]})) {
						$cpdhits->{rep}->{$cpds->[$i]->id()}->{$peakid} = 1;
						$cpdhits->{all}->{$cpds->[$i]->id()}->{$peakid} = 1;
					}
				}
			}
			$combined->add("modelcompounds",$cpds->[$i]->cloneObject());
		} 
	}
	if (defined($original_cpd_hash->{$cpds->[$i]->inchikey()}->{damage})) {
		$cpds->[$i]->dblinks()->{olddmg} = $original_cpd_hash->{$cpds->[$i]->inchikey()}->{damage};
		$original_cpd_id_hash->{damage}->{$rxns->[$i]->id()} = 1;
	}
	if (defined($original_cpd_hash->{$cpds->[$i]->inchikey()}->{repair})) {
		$cpds->[$i]->dblinks()->{oldrep} = $original_cpd_hash->{$cpds->[$i]->inchikey()}->{repair};
		$original_cpd_id_hash->{repair}->{$rxns->[$i]->id()} = 1;
	}
}

print "Total repair compounds:".$repcount."\n";
print "Repair peaks:".keys(%{$peakhits->{rep}})."\n";
print "Repair seed compounds:".keys(%{$cpdhits->{rep}})."\n";
print "Repair other compounds:".keys(%{$otherhits->{rep}})."\n";
print "Total peaks:".keys(%{$peakhits->{all}})."\n";
print "Total seed compounds:".keys(%{$cpdhits->{all}})."\n";
print "Total other compounds:".keys(%{$otherhits->{all}})."\n";

foreach my $cpdid (keys(%{$cpdhits->{db}})) {
	foreach my $peakid (keys(%{$cpdhits->{db}->{$cpdid}})) {
		$cpdhits->{all}->{$cpdid}->{$peakid} = 1;
		$peakhits->{all}->{$peakid}->{$cpdid} = 1;
		$peakhits->{db}->{$peakid}->{$cpdid} = 1;
	}
}
print "DB peaks:".keys(%{$peakhits->{db}})."\n";
print "DB seed compounds:".keys(%{$cpdhits->{db}})."\n";
print "Total peaks:".keys(%{$peakhits->{all}})."\n";
print "Total seed compounds:".keys(%{$cpdhits->{all}})."\n";
print "Total other compounds:".keys(%{$otherhits->{all}})."\n";

$cpds = $combined->modelcompounds();
my $combined_hash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	$combined_hash->{$cpds->[$i]->id()} = 1;
}

$repcount = 0;
for (my $i=0; $i < @{$rxns}; $i++) {
	my $id = $rxns->[$i]->id();	
	if ($id =~ m/pkr/) {
		$id =~ s/pkr/rep/;
		$rxns->[$i]->id($id);
	}
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage})) {
		$rxns->[$i]->dblinks()->{olddmg} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{damage};
		$original_rxn_id_hash->{damage}->{$rxns->[$i]->id()} = 1;
	}
	if (defined($original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair})) {
		$rxns->[$i]->dblinks()->{oldrep} = $original_rxn_hash->{$rxns->[$i]->createEquation({indecies => 0,format=>"inchikey",hashed=>1,protons=>0,direction=>0})}->{repair};
		$original_rxn_id_hash->{repair}->{$rxns->[$i]->id()} = 1;
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
		$repcount++;
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

print "Total repair reactions:".$repcount."\n";

#my $wsmeta = $impl->util_save_object($dam_model,"29280/MMSyn3_damage.fix");
#$wsmeta = $impl->util_save_object($rep_model,"29280/MMSyn3_repair.fix");

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
my $wsmeta = $impl->util_save_object($combined,"29280/MMSyn3_combined.fix");

#Now adding original model to combined model
$cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (!defined($combined->getObject("modelcompounds",$cpds->[$i]->id()))) {
		$combined->add("modelcompounds",$cpds->[$i]->cloneObject());
	}
}
$rxns = $model->modelreactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$combined->add("modelreactions",$rxns->[$i]->cloneObject());
}

#Saving fully combined model
my $wsmeta = $impl->util_save_object($combined,"29280/FullCombined.fix");

my $match_dmg_cpd = [0,0];
my $match_rep_cpd = [0,0];
my $match_dmg_rxn = [0,0];
my $match_rep_rxn = [0,0];
foreach my $cpd (keys(%{$original_cpd_id_hash->{damage}})) {
	if ($original_cpd_id_hash->{damage}->{$cpd} == 0) {
		$match_dmg_cpd->[0]++;
	} else {
		$match_dmg_cpd->[1]++;
	}
}
foreach my $cpd (keys(%{$original_cpd_id_hash->{repair}})) {	
	if ($original_cpd_id_hash->{repair}->{$cpd} == 0) {
		$match_rep_cpd->[0]++;
	} else {
		$match_rep_cpd->[1]++;
	}
}
foreach my $rxn (keys(%{$original_rxn_id_hash->{damage}})) {
	if ($original_rxn_id_hash->{damage}->{$rxn} == 0) {
		$match_dmg_rxn->[0]++;
	} else {
		$match_dmg_rxn->[1]++;
	}
}
foreach my $rxn (keys(%{$original_rxn_id_hash->{repair}})) {
	if ($original_rxn_id_hash->{repair}->{$rxn} == 0) {
		$match_dmg_rxn->[0]++;
	} else {
		$match_dmg_rxn->[1]++;
	}
}

print "Repair compounds:".$match_rep_cpd->[0]."/".$match_rep_cpd->[1]."\n";
print "Damage compounds:".$match_dmg_cpd->[0]."/".$match_dmg_cpd->[1]."\n";
print "Repair reactions:".$match_rep_rxn->[0]."/".$match_rep_rxn->[1]."\n";
print "Damage reactions:".$match_dmg_rxn->[0]."/".$match_dmg_rxn->[1]."\n";