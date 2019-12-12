use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("CombinedModel",$workspace));
my $cpds = $model->modelcompounds();

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/Peaklist.txt");
my $peakstring = "";
my $peakids;
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	#First check for ModelSEED ID, then check for inchi-key
	if (defined($model->getObject("modelcompounds",$array->[1]))) {
		if (length($peakstring) > 0) {
			$peakstring .= ";";
		}
		$peakstring .= "peak.".$i.":".$array->[1];
		$peakids->{$array->[1]} = 1;
	} else {
		#Then check if there is an inchi string match
		for (my $j=0; $j < @{$cpds}; $j++) {
			if ($cpds->[$j]->inchikey() eq $array->[0]) {
				if (length($peakstring) > 0) {
					$peakstring .= ";";
				}
				$peakstring .= "peak.".$i.":".$cpds->[$j]->id();
				$peakids->{$cpds->[$j]->id()} = 1;
				last;
			}
		}
	}
}

#print $peakstring."\n";

#my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
#	fbamodel_id => "iMB155.trans.noex.fix",
#	workspace => $workspace,
#	media_id => "Complete",
#	media_workspace => "KBaseMedia",
#	source_fbamodel_id => "CombinedModel",
#	target_reaction => "bio1",
#	fbamodel_output_id => "iMB155.trans.noex.fix.metgf",
#	metabolomics_peak_data => $peakstring
#});

exit;

my $gfmodel = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.trans.noex.fix.metgf",$workspace));

my $rxnhash;
my $peaks;
my $rxns = $gfmodel->modelreactions();
my $gfrxns = [];
my $gfrxn_hash = {};
for (my $i=0; $i < @{$rxns}; $i++) {
	$peaks->{$rxns->[$i]->id()} = 0;
	$rxnhash->{$rxns->[$i]->id()} = $rxns->[$i];
	my $gfdata = $rxns->[$i]->gapfill_data();
	my $gfnum = keys(%{$gfdata});
	if ($gfnum > 0) {
		$gfrxn_hash->{$rxns->[$i]->id()} = 1;
		push(@{$gfrxns},$rxns->[$i]);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $k=0; $k < @{$rgts}; $k++) {
		if (defined($peakids->{$rgts->[$k]->modelcompound()->id()})) {
			$peaks->{$rxns->[$i]->id()}++;
		}
	}
}

my $marks = {};
my $rxn_marks = {};
my $cofactor_hash;
my $cofactors = Bio::KBase::constants::cofactors();
for (my $i=0; $i < @{$cofactors}; $i++) {
	$cofactor_hash->{$cofactors->[$i]} = 1;
}

my $sets;
my $setid = 0;
print "SetID\tRound\tModelRxn\tID\tEquation\tDefinition\tPeaks\tNewPeakScore\tModelPeakScore\tTotalPeakScore\n";
for (my $i=0; $i < @{$gfrxns}; $i++) {
	if (!defined($rxn_marks->{$gfrxns->[$i]->id()})) {
		$setid++;
		my $peakcount = $peaks->{$gfrxns->[$i]->id()};
		my $modelpeakcount = 0;
		$rxn_marks->{$gfrxns->[$i]->id()} = 1;
		my $round = 0;
		my $current_set = {$gfrxns->[$i]->id() => 0};
		my $current_model_rxn = {};
		my $gfrxn = $gfrxns->[$i];
		my $rgts = $gfrxn->modelReactionReagents();
		for (my $j=0; $j < @{$rgts}; $j++) {
			if ($rgts->[$j]->modelcompound()->id() =~ m/(cpd\d+)/) {
				my $msid = $1;
				if (!defined($cofactor_hash->{$msid})) {
					$marks->{$rgts->[$j]->modelcompound()->id()} = 1;
				}
			}
		}
		my $new = 1;
		while ($new == 1) {
			$new = 0;
			$round++;
			for (my $j=0; $j < @{$rxns}; $j++) {
				if (!defined($rxn_marks->{$rxns->[$j]->id()})) {
					#Checking if reaction contains a marked reagent
					my $rgts = $rxns->[$j]->modelReactionReagents();
					my $found = 0;
					for (my $k=0; $k < @{$rgts}; $k++) {
						if (defined($marks->{$rgts->[$k]->modelcompound()->id()})) {
							$found = 1;
							last;
						}
					}
					if ($found == 1) {
						if (defined($gfrxn_hash->{$rxns->[$j]->id()})) {
							$current_set->{$rxns->[$j]->id()} = $round;
							$peakcount += $peaks->{$rxns->[$j]->id()};
							for (my $k=0; $k < @{$rgts}; $k++) {
								if ($rgts->[$k]->modelcompound()->id() =~ m/(cpd\d+)/) {
									my $msid = $1;
									if (!defined($cofactor_hash->{$msid})) {
										$marks->{$rgts->[$k]->modelcompound()->id()} = 1;
									}
								}
							}
							$new = 1;
						} else {
							$modelpeakcount += $peaks->{$rxns->[$j]->id()};
							$current_model_rxn->{$rxns->[$j]->id()} = $round;
						}
						$rxn_marks->{$rxns->[$j]->id()} = 1;
					}
				}
			}
		}
		#Unmarking model reactions and their reagents
		for (my $j=0; $j < @{$rxns}; $j++) {
			if (!defined($gfrxn_hash->{$rxns->[$j]->id()})) {
				delete $rxn_marks->{$rxns->[$j]->id()};
				my $rgts = $rxns->[$j]->modelReactionReagents();
				for (my $k=0; $k < @{$rgts}; $k++) {
					if ($rgts->[$k]->modelcompound()->id() =~ m/(cpd\d+)/) {
						my $msid = $1;
						if (!defined($cofactor_hash->{$msid})) {
							delete $marks->{$rgts->[$k]->modelcompound()->id()};
						}
					}
				}
			}
		}
		#Printing reaction set
		foreach my $rxnid (keys(%{$current_set})) {
			print $setid."\t".$current_set->{$rxnid}."\tPredicted\t".$rxnid."\t".$rxnhash->{$rxnid}->equation()."\t".$rxnhash->{$rxnid}->definition()."\t".$peaks->{$rxnid}."\t".$peakcount."\t".$modelpeakcount."\t".($peakcount+$modelpeakcount)."\n";
		}
		foreach my $rxnid (keys(%{$current_model_rxn})) {
			print $setid."\t".$current_model_rxn->{$rxnid}."\tModel\t".$rxnid."\t".$rxnhash->{$rxnid}->equation()."\t".$rxnhash->{$rxnid}->definition()."\t".$peaks->{$rxnid}."\t".$peakcount."\t".$modelpeakcount."\t".($peakcount+$modelpeakcount)."\n";
		}
	}
}

