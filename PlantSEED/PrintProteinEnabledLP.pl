use strict;
use Data::Dumper;
use JSON::XS;
use JSON;

local $| = 1;

my $modelfile = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ZMayz1.json";
#my $modelfile = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ZMayz2.json";
#my $modelfile = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ZMayz16.json";
my $filename = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MaizeFBA1.lp";
my $kcat_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionKCatFlux.txt";
my $protein_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionProtein.txt";
my $measured_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MeasuredReaction.txt";
my $compound_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/Compounds.json";
#Run just for leaf on it's own
my $datasets = {Mature_leaf => 0};
#Run for leaf and root
#my $datasets = {Mature_leaf => 1,root_MZ=>2};
#my $datasets = {Mature_leaf => 1,root_MZ=>2};
my $protein_fraction = {
	root_MZ => 0,
	root_EZ => 0,
	root_Cortex => 0,
	root_Stele => 0,
	Embryo_20_DAP => 0,
	Embryo_38_DAP => 0,
	Endosperm_8_DAP => 0,
	Endosperm_10_DAP => 0,
	Endosperm_12_DAP => 0,
	Endosperm_Crown_27_DAP => 0,
	Pericarp_Aleurone_27_DAP => 0,
	GermEmbryo_2_DAI => 0,
	Zone_1 => 0,
	Zone_2 => 0,
	Zone_3 => 0,
	Mature_leaf => 0.0948
};

my $vopt_coef = 1;
my $measured_coef = 100;
my $kfit_coef = 1;
my $excluded_coef = 1;

#Loading the model
my $data = [];
open (my $fh, "<", $modelfile) || print "Couldn't open $modelfile: $!";
while (my $Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	push(@{$data},$Line);
}
close($fh);

#Converting model file to data structure
my $model = decode_json join("\n",@{$data});

#Loading vopt
my $reaction_vopt_hash = {};
open ($fh, "<", $kcat_file) || print "Couldn't open $kcat_file: $!";
my $heading;
while (my $Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	if (!defined($heading)) {
		$heading = [split(/\t/,$Line)];
	} else {
		my $array = [split(/\t/,$Line)];
		for (my $i=1; $i < @{$array}; $i++) {
			$reaction_vopt_hash->{$array->[0]}->{$heading->[$i]} = $array->[$i];
		}
	}
}
close($fh);

#Loading reaction proteins
my $reaction_protein_hash = {};
open ($fh, "<", $protein_file) || print "Couldn't open $protein_file: $!";
$heading = undef;
while (my $Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	if (!defined($heading)) {
		$heading = [split(/\t/,$Line)];
	} else {
		my $array = [split(/\t/,$Line)];
		for (my $i=1; $i < @{$array}; $i++) {
			$reaction_protein_hash->{$array->[0]}->{$heading->[$i]} = $array->[$i];
		}
	}
}
close($fh);

#Loading reaction measured flux
my $reaction_measured_hash = {};
open ($fh, "<", $measured_file) || print "Couldn't open $measured_file: $!";
my $Line = <$fh>;
while ($Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	my $array = [split(/\t/,$Line)];
	$reaction_measured_hash->{$array->[0]}->{$array->[1]} = $array->[2];
}
close($fh);

#Loading compound hash with molecular weight information
my $compound_hash = {};
$data = "";
open ($fh, "<", $compound_file) || print "Couldn't open $compound_file: $!";
while ($Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	$data .= $Line."\n";
}
my $cpddata = decode_json $data;
for (my $i=0; $i < @{$cpddata}; $i++) {
	$compound_hash->{$cpddata->[$i]->{id}} = $cpddata->[$i];
}

#Initializing an empty problem
my $problem = {
	varibles => [],
	constraints => [],
	objective => {
		variables => [],
		coefficients => [],
		quadratic => []
	}
};

my $metabolites = {};
my $bio_metabolites = {};
my $bios = $model->{biomasses};
my $biohash = {};
foreach my $bio (@{$bios}) {
	$biohash->{$bio->{id}} = $bio;
}
foreach my $dataset (keys(%{$datasets})) {
	my $index = $datasets->{$dataset};
	my $bio = $biohash->{"bio".$index};
	my $rgts = $bio->{biomasscompounds};
	foreach my $rgt (@{$rgts}) {
		if ($rgt->{coefficient} < 0) {
			my $cpdid = $rgt->{modelcompound_ref};
			$cpdid =~ s/.+\///g;
			$bio_metabolites->{$bio->{id}}->{$cpdid} = 1;
			$metabolites->{$cpdid} = $rgt->{coefficient};
		}
	}
}

#Creating compound variables and constraints
my $cpd_obj = {};
my $cpds = $model->{modelcompounds};
foreach my $cpd (@{$cpds}) {
	#Creating mass balanced constraints for every single compound
	$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance} = {
		name => "MB_".$cpd->{id},
		type => "MassBalance",
		variables => [],
		coefficients => [],
		rhs => 0,
		sign => "="
	};
	push(@{$problem->{constraints}},$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance});
	#Creating a drain flux for every extracellular compound - drain fluxes look like this: => cpdXXXXX 
	if ($cpd->{id} =~ m/_e\d+$/ || $cpd->{id} =~ m/cpd11416/ || $cpd->{id} =~ m/cpd02701/ ||  defined($metabolites->{$cpd->{id}})) {
		if (defined($metabolites->{$cpd->{id}})) { #if biomass drain flux
			$cpd_obj->{$cpd->{id}}->{variables}->{Flux} = {
				type => "Flux",
				name => "F_drain_".$cpd->{id},
				binary => 0,
				upperbound => 1000,
				lowerbound => -1000
			};
		} else {
			#If we want to simulate something other than complete media, we need to change this
			$cpd_obj->{$cpd->{id}}->{variables}->{Flux} = {
				type => "Flux",
				name => "F_".$cpd->{id},
				binary => 0,
				upperbound => 1000,
				lowerbound => -1000
			};
		}
		if ($cpd->{id} =~ m/cpd11416/ || $cpd->{id} =~ m/cpd02701/) {
			$cpd_obj->{$cpd->{id}}->{variables}->{Flux}->{upperbound} = 0;
		}
		push(@{$problem->{variables}},$cpd_obj->{$cpd->{id}}->{variables}->{Flux});
		#Adding drain flux to mass balanced constraint
		push(@{$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance}->{variables}},$cpd_obj->{$cpd->{id}}->{variables}->{Flux});
		push(@{$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance}->{coefficients}},1);

	}
}

#Creating compound variables and constraints
my $rxn_obj = {};
my $rxns = $model->{modelreactions};
foreach my $rxn (@{$rxns}) {
	#Creating flux variables for every reaction with bounds set according to directionality
	$rxn_obj->{$rxn->{id}}->{variables}->{Flux} = {
		type => "Flux",
		name => "F_".$rxn->{id},
		binary => 0,
		upperbound => 1000,
		lowerbound => -1000
	};
	if ($rxn->{direction} eq ">") {
		$rxn_obj->{$rxn->{id}}->{variables}->{Flux}->{lowerbound} = 0;
	} elsif ($rxn->{direction} eq "<") {
		$rxn_obj->{$rxn->{id}}->{variables}->{Flux}->{upperbound} = 0;
	}
	push(@{$problem->{variables}},$rxn_obj->{$rxn->{id}}->{variables}->{Flux});
	#Adding reaction flux to associated mass balance constraints
	my $rgts = $rxn->{modelReactionReagents};
	foreach my $rgt (@{$rgts}) {
		my $cpdid = $rgt->{modelcompound_ref};
		$cpdid =~ s/.+\///g;
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{variables}},$rxn_obj->{$rxn->{id}}->{variables}->{Flux});
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{coefficients}},$rgt->{coefficient});
	}
}

#Creating biomass variables and constraints
my $bios = $model->{biomasses};
foreach my $bio (@{$bios}) {
	#Creating flux variables for every biomass reaction
	$rxn_obj->{$bio->{id}}->{variables}->{Flux} = {
		type => "Flux",
		name => "F_".$bio->{id},
		binary => 0,
		upperbound => 1000,
		lowerbound => 0
	};
	push(@{$problem->{variables}},$rxn_obj->{$bio->{id}}->{variables}->{Flux});
	#Adding reaction flux to associated mass balance constraints
	my $rgts = $bio->{biomasscompounds};
	foreach my $rgt (@{$rgts}) {
		my $cpdid = $rgt->{modelcompound_ref};
		$cpdid =~ s/.+\///g;
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{variables}},$rxn_obj->{$bio->{id}}->{variables}->{Flux});
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{coefficients}},$rgt->{coefficient});
	}
}

#Adding flexible biomass drain constraints
foreach my $dataset (keys(%{$datasets})) {
	my $index = $datasets->{$dataset};
	$rxn_obj->{"bio".$index}->{constraints}->{FlexibleBiomass} = {
		name => "Flex_bio".$index,
		type => "FlexibleBiomassConstraint",
		variables => [],
		coefficients => [],
		rhs => 0,
		sign => "="
	};
	if (defined($bio_metabolites->{"bio".$index})) {
		foreach my $met (keys(%{$bio_metabolites->{"bio".$index}})) {
			my $mass = 1;
			if ($met =~ m/(cpd\d+)/ && defined($compound_hash->{$1})) {
				$mass = $compound_hash->{$1}->{mass};
			}
			push(@{$rxn_obj->{"bio".$index}->{constraints}->{FlexibleBiomass}->{variables}},$cpd_obj->{$met}->{variables}->{Flux});
			push(@{$rxn_obj->{"bio".$index}->{constraints}->{FlexibleBiomass}->{coefficients}},$mass);
			push(@{$problem->{constraints}},{
				name => "BCMax_".$met,
				type => "BiomassCompoundConstraint",
				variables => [$rxn_obj->{"bio".$index}->{variables}->{Flux},$cpd_obj->{$met}->{variables}->{Flux}],
				coefficients => [0.75 * abs($metabolites->{$met}),-1],
				rhs => 0,
				sign => ">="
			});
			push(@{$problem->{constraints}},{
				name => "BCMin_".$met,
				type => "BiomassCompoundConstraint",
				variables => [$rxn_obj->{"bio".$index}->{variables}->{Flux},$cpd_obj->{$met}->{variables}->{Flux}],
				coefficients => [-0.75 * abs($metabolites->{$met}),-1],
				rhs => 0,
				sign => "<="
			});
		}
	}
}

#Creating objective
my $rxnhash = {};
#Adding measured flux constraint
foreach my $rxnid (keys(%{$reaction_measured_hash})) {
	#Adding constraint: v - vm = vfit
	foreach my $dataset (keys(%{$datasets})) {
		my $full_id = $rxnid.$datasets->{$dataset};
		$rxnhash->{$full_id} = 1;
		if (defined($rxn_obj->{$full_id}->{variables}->{Flux}) || $cpd_obj->{$full_id}->{variables}->{Flux}) {
			my $fluxvar;
			if (defined($rxn_obj->{$full_id}->{variables}->{Flux})) {
				$fluxvar = $rxn_obj->{$full_id}->{variables}->{Flux};
			} else {
				$fluxvar = $cpd_obj->{$full_id}->{variables}->{Flux};
			}
			#Adding fit variable
			$rxn_obj->{$full_id}->{variables}->{Measured} = {
				type => "Measured",
				name => "M_".$full_id,
				binary => 0,
				upperbound => 1000,
				lowerbound => -1000
			};
			push(@{$problem->{variables}},$rxn_obj->{$full_id}->{variables}->{Measured});
			$rxn_obj->{$full_id}->{constraints}->{Measured} = {
				name => "MC_".$full_id,
				type => "Measured",
				variables => [$fluxvar,$rxn_obj->{$full_id}->{variables}->{Measured}],
				coefficients => [1,-1],
				rhs => $reaction_measured_hash->{$rxnid}->{$dataset},
				sign => "="
			};
			push(@{$problem->{constraints}},$rxn_obj->{$full_id}->{constraints}->{Measured});
			push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{Measured});
			push(@{$problem->{objective}->{coefficients}},$measured_coef);
			push(@{$problem->{objective}->{quadratic}},1);
		}
	}
}

#Adding kcat optimal flux constraint
foreach my $rxnid (keys(%{$reaction_vopt_hash})) {
	if (!defined($rxnhash->{$rxnid})) {
		#Adding constraint: v - vm = vfit
		foreach my $dataset (keys(%{$datasets})) {
			my $full_id = $rxnid.$datasets->{$dataset};
			$rxnhash->{$full_id} = 1;
			if (defined($rxn_obj->{$full_id}->{variables}->{Flux})) {
				#Adding fit variable
				$rxn_obj->{$full_id}->{variables}->{KcatOpt} = {
					type => "KcatOpt",
					name => "K_".$full_id,
					binary => 0,
					upperbound => 1000,
					lowerbound => -1000
				};
				push(@{$problem->{variables}},$rxn_obj->{$full_id}->{variables}->{KcatOpt});
				$rxn_obj->{$full_id}->{constraints}->{KcatOpt} = {
					name => "KC_".$full_id,
					type => "KcatOpt",
					variables => [$rxn_obj->{$full_id}->{variables}->{Flux},$rxn_obj->{$full_id}->{variables}->{KcatOpt}],
					coefficients => [1,-1],
					rhs => $protein_fraction->{$dataset}*$reaction_vopt_hash->{$rxnid}->{$dataset},
					sign => "="
				};
				push(@{$problem->{constraints}},$rxn_obj->{$full_id}->{constraints}->{KcatOpt});
				push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{KcatOpt});
				push(@{$problem->{objective}->{coefficients}},$vopt_coef);
				push(@{$problem->{objective}->{quadratic}},1);
			}
		}
	}
}

#Adding kfit optimal flux constraint
my $numdatasets = keys(%{$datasets});
foreach my $rxnid (keys(%{$reaction_protein_hash})) {
	if (!defined($rxnhash->{$rxnid})) {
		#Adding constraint: v - vm = vfit
		foreach my $dataset (keys(%{$datasets})) {
			my $full_id = $rxnid.$datasets->{$dataset};
			$rxnhash->{$full_id} = 1;
			#Adding kfit variable
			if (defined($rxn_obj->{$full_id}->{variables}->{Flux})) {
				if (!defined($reaction_protein_hash->{$rxnid}->{$dataset}) || $reaction_protein_hash->{$rxnid}->{$dataset} == 0 || $numdatasets == 1) {
					push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{Flux});
				} else {
					$rxn_obj->{$full_id}->{variables}->{KValue} = {
						type => "KValue",
						name => "KV_".$full_id,
						binary => 0,
						upperbound => 1000000,
						lowerbound => -1000000
					};
					push(@{$problem->{variables}},$rxn_obj->{$full_id}->{variables}->{KValue});
					#Adding fit variable
					$rxn_obj->{$full_id}->{variables}->{KfitOpt} = {
						type => "KfitOpt",
						name => "KF_".$full_id,
						binary => 0,
						upperbound => 1000,
						lowerbound => -1000
					};
					push(@{$problem->{variables}},$rxn_obj->{$full_id}->{variables}->{KfitOpt});
					$rxn_obj->{$full_id}->{constraints}->{KfitOpt} = {
						name => "KF_".$full_id,
						type => "KfitOpt",
						variables => [$rxn_obj->{$full_id}->{variables}->{Flux},$rxn_obj->{$full_id}->{variables}->{KfitOpt},$rxn_obj->{$full_id}->{variables}->{KValue}],
						coefficients => [1,-1,-1*$protein_fraction->{$dataset}*$reaction_protein_hash->{$rxnid}->{$dataset}],
						rhs => 0,
						sign => "="
					};
					push(@{$problem->{constraints}},$rxn_obj->{$full_id}->{constraints}->{KfitOpt});
					push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{KfitOpt});
				}
				push(@{$problem->{objective}->{coefficients}},$vopt_coef);
				push(@{$problem->{objective}->{quadratic}},1);
			}
		}
	}
}

#Adding all reactions with excluded flux to objective to be set to zero
foreach my $rxn (@{$rxns}) {
	if (!defined($rxnhash->{$rxn->{id}})) {
		push(@{$problem->{objective}->{variables}},$rxn_obj->{$rxn->{id}}->{variables}->{Flux});
		push(@{$problem->{objective}->{coefficients}},$vopt_coef);
		push(@{$problem->{objective}->{quadratic}},1);
	}	
}

##Replacing objective with biomass
#$problem->{objective}->{variables} = [$rxn_obj->{bio0}->{variables}->{Flux}];
#$problem->{objective}->{coefficients} = [-1];
#$problem->{objective}->{quadratic} = [0];

#Printing LP file
my $output = ['\* Problem: ProteomDrivenModelFBA *\\',"","Minimize","","","Subject To","","","Bounds","","","Binary","","","End"];
#Printing objective
$output->[3] .= "obj:";
#Printing nonquadratic terms first
my $count = 0;
for (my $i=0; $i < @{$problem->{objective}->{variables}}; $i++) {
	if ($problem->{objective}->{quadratic}->[$i] == 0) {
		$count++;
		if ($count % 6 == 5) {
			$output->[3] .= "\n";
		}
		if ($problem->{objective}->{coefficients}->[$i] > 0) {
			$output->[3] .= " + ".$problem->{objective}->{coefficients}->[$i]." ".$problem->{objective}->{variables}->[$i]->{name};
		} elsif ($problem->{objective}->{coefficients}->[$i] < 0) {
			$output->[3] .= " - ".(-1*$problem->{objective}->{coefficients}->[$i])." ".$problem->{objective}->{variables}->[$i]->{name};
		}
	}
}
#Now printing quadratic terms
$count = 0;
for (my $i=0; $i < @{$problem->{objective}->{variables}}; $i++) {
	if ($problem->{objective}->{quadratic}->[$i] == 1) {
		if ($count == 0) {
			if ($count > 0) {
				$output->[3] .= " +";
			}
			$output->[3] .= " ["; 
		}
		$count++;
		if ($count % 6 == 5) {
			$output->[3] .= "\n";
		}
		if ($problem->{objective}->{coefficients}->[$i] > 0) {
			$output->[3] .= " + ".$problem->{objective}->{coefficients}->[$i]." ".$problem->{objective}->{variables}->[$i]->{name}."^2";
		} elsif ($problem->{objective}->{coefficients}->[$i] < 0) {
			$output->[3] .= " - ".(-1*$problem->{objective}->{coefficients}->[$i])." ".$problem->{objective}->{variables}->[$i]->{name}."^2";
		}
	}
}
if ($count > 0) {
	$output->[3] .= "]/2";
}
#Printing constraints
for (my $j=0; $j < @{$problem->{constraints}}; $j++) {
	if ($j > 0) {
		$output->[6] .= "\n";
	}
	my $constraint = $problem->{constraints}->[$j];
	$output->[6] .= " ".$constraint->{name}.":";
	for (my $i=0; $i < @{$constraint->{variables}}; $i++) {
		if ($i % 6 == 5) {
			$output->[6] .= "\n";
		}
		if ($constraint->{coefficients}->[$i] > 0) {
			$output->[6] .= " + ".$constraint->{coefficients}->[$i]." ".$constraint->{variables}->[$i]->{name};
		} elsif ($constraint->{coefficients}->[$i] < 0) {
			$output->[6] .= " - ".(-1*$constraint->{coefficients}->[$i])." ".$constraint->{variables}->[$i]->{name};
		}
	}
	$output->[6] .= " ".$constraint->{sign}." ".$constraint->{rhs};	
}

#Printing bounds
for (my $i=0; $i < @{$problem->{variables}}; $i++) {
	$output->[9] .= $problem->{variables}->[$i]->{lowerbound}." <= ".$problem->{variables}->[$i]->{name}." <= ".$problem->{variables}->[$i]->{upperbound}."\n";
}
#Printing binary list
my $count = 0;
for (my $i=0; $i < @{$problem->{variables}}; $i++) {
	if ($problem->{variables}->[$i]->{binary} == 1) {
		if ($count % 8 == 7) {
			$output->[12] .= "\n";
		}
		$count++;
		$output->[12] .= " ".$problem->{variables}->[$i]->{name};
	}
}
#Printing file
open ( my $fout, ">", $filename) || print "Failure to open file: $filename, $!";
foreach my $Item (@{$output}) {
	print $fout $Item."\n";
}
close($fout);