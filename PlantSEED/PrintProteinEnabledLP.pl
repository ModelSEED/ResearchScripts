use strict;
use Data::Dumper;
use JSON::XS;
use JSON;

local $| = 1;

my $modelfile = $ARGV[0];
my $filename = $ARGV[1];
$modelfile = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ZMayz.json";
#Converting model file to data structure
my $model = decode_json join("\n",@{$data});

$filename = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MaizeFBA.lp";
my $kcat_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionKCatFlux.txt";
my $protein_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionProtein.txt";
my $measured_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MeasuredReaction.txt";
my $datasets = {};

my $vopt_coef = 1;
my $measured_coef = 1;
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



#Initializing an empty problem
my $problem = {
	varibles => [],
	constraints => [],
	objective => {
		variables => [],
		coefficients => [-1],
		quadratic => [0]
	}
};

#TODO : Add the variables and constraints for flexible biomass
#TODO : Add the additional variables and constraints for the kinetic formulation

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
	if ($cpd->{id} =~ m/_e\d+$/ || $cpd->{id} =~ m/cpd11416/ || $cpd->{id} =~ m/cpd02701/) {
		#If we want to simulate something other than complete media, we need to change this
		$cpd_obj->{$cpd->{id}}->{variables}->{Flux} = {
			type => "Flux",
			name => "F_".$cpd->{id},
			binary => 0,
			upperbound => 1000,
			lowerbound => -1000
		};
		if ($cpd->{id} =~ m/cpd11416/ || $cpd->{id} =~ m/cpd02701/) {
			$cpd_obj->{$cpd->{id}}->{variables}->{Flux}->{upperbound} = 0;
		}
		push(@{$problem->{variables}},$cpd_obj->{$cpd->{id}}->{variables}->{Flux});
		#Adding drain flux to mass balanced constraint
		push(@{$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance}->{variables}},$cpd_obj->{$cpd->{id}}->{variables}->{Flux});
		push(@{$cpd_obj->{$cpd->{id}}->{constraints}->{MassBalance}->{coefficients}},1);
	}
}

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

my $bio_obj = {};
my $bios = $model->{biomasses};
foreach my $bio (@{$bios}) {
	#Creating flux variables for every biomass reaction
	$bio_obj->{$bio->{id}}->{variables}->{Flux} = {
		type => "Flux",
		name => "F_".$bio->{id},
		binary => 0,
		upperbound => 1000,
		lowerbound => 0
	};
	push(@{$problem->{variables}},$bio_obj->{$bio->{id}}->{variables}->{Flux});
	#Adding reaction flux to associated mass balance constraints
	my $rgts = $bio->{biomasscompounds};
	foreach my $rgt (@{$rgts}) {
		my $cpdid = $rgt->{modelcompound_ref};
		$cpdid =~ s/.+\///g;
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{variables}},$bio_obj->{$bio->{id}}->{variables}->{Flux});
		push(@{$cpd_obj->{$cpdid}->{constraints}->{MassBalance}->{coefficients}},$rgt->{coefficient});
	}
}

#Creating objective
my $rxnhash = {};
#Adding measured flux constraint
foreach my $rxnid (keys(%{$reaction_measured_hash})) {
	#Adding constraint: v - vm = vfit
	foreach my $dataset (keys(%{$reaction_measured_hash->{$rxnid}})) {
		my $full_id = $rxnid.$datasets->{$dataset};
		$rxnhash->{$full_id} = 1;
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
			variables => [$rxn_obj->{$full_id}->{variables}->{Flux},$rxn_obj->{$full_id}->{variables}->{Measured}],
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
#Adding kcat optimal flux constraint
foreach my $rxnid (keys(%{$reaction_vopt_hash})) {
	if (!defined($rxnhash->{$rxnid})) {
		#Adding constraint: v - vm = vfit
		foreach my $dataset (keys(%{$reaction_vopt_hash->{$rxnid}})) {
			my $full_id = $rxnid.$datasets->{$dataset};
			$rxnhash->{$full_id} = 1;
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
				rhs => $reaction_vopt_hash->{$rxnid}->{$dataset},
				sign => "="
			};
			push(@{$problem->{constraints}},$rxn_obj->{$full_id}->{constraints}->{KcatOpt});
			push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{KcatOpt});
			push(@{$problem->{objective}->{coefficients}},$vopt_coef);
			push(@{$problem->{objective}->{quadratic}},1);
		}
	}
}
#Adding kfit optimal flux constraint
foreach my $rxnid (keys(%{$reaction_protein_hash})) {
	if (!defined($rxnhash->{$rxnid})) {
		#Adding constraint: v - vm = vfit
		foreach my $dataset (keys(%{$reaction_protein_hash->{$rxnid}})) {
			my $full_id = $rxnid.$datasets->{$dataset};
			$rxnhash->{$full_id} = 1;
			#Adding kfit variable
			$rxn_obj->{$full_id}->{variables}->{KValue} = {
				type => "KValue",
				name => "KV_".$full_id,
				binary => 0,
				upperbound => 1000,
				lowerbound => -1000
			};
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
				variables => [$rxn_obj->{$full_id}->{variables}->{Flux},$rxn_obj->{$full_id}->{variables}->{KfitOpt}],
				coefficients => [1,-1,-1*$reaction_protein_hash->{$rxnid}->{$dataset}],
				rhs => 0,
				sign => "="
			};
			push(@{$problem->{constraints}},$rxn_obj->{$full_id}->{constraints}->{KfitOpt});
			push(@{$problem->{objective}->{variables}},$rxn_obj->{$full_id}->{variables}->{KfitOpt});
			push(@{$problem->{objective}->{coefficients}},$vopt_coef);
			push(@{$problem->{objective}->{quadratic}},1);
		}
	}
}

my $kfit_coef = 1;
my $excluded_coef = 1;

$problem->{objective}->{variables} = [$bio_obj->{"bio1"}->{variables}->{Flux}];

#Printing LP file
my $output = ['\* Problem: ProteomDrivenModelFBA *\\',"","Minimize","","","Subject To","","","Bounds","","","Binary","","","End"];
#Printing objective
$output->[3] .= "obj:";
for (my $i=0; $i < @{$problem->{objective}->{variables}}; $i++) {
	if ($i % 6 == 5) {
		$output->[3] .= "\n";
	}
	if ($problem->{objective}->{coefficients}->[$i] > 0) {
		$output->[3] .= " + ".$problem->{objective}->{coefficients}->[$i]." ".$problem->{objective}->{variables}->[$i]->{name};
	} elsif ($problem->{objective}->{coefficients}->[$i] < 0) {
		$output->[3] .= " - ".(-1*$problem->{objective}->{coefficients}->[$i])." ".$problem->{objective}->{variables}->[$i]->{name};
	}
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