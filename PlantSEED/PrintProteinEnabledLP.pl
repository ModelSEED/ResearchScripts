use strict;
use Data::Dumper;
use JSON::XS;
use JSON;

local $| = 1;

my $modelfile = $ARGV[0];
my $filename = $ARGV[1];
$modelfile = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ZMayz.json";
$filename = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MaizeFBA.lp";
my $kcat_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionKCatFlux.txt";
my $protein_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/ReactionProtein.txt";
my $measured_file = "/Users/chenry/Dropbox/workspace/PlantSEED/ProteomeModeling/MeasuredReaction.txt";
my $datasets = {Mature_leaf => 1,root_MZ=>2};
my $protein_fraction = {
	
};

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

#Loading reaction proteins
my $reaction_measured_hash = {};
open ($fh, "<", $measured_file) || print "Couldn't open $measured_file: $!";
while (my $Line = <$fh>) {
	$Line =~ s/\r//;
	chomp($Line);
	my $array = [split(/\t/,$Line)];
	$reaction_measured_hash->{$array->[0]}->{$array->[1]} = $array->[2];
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

#my @metabolites=('cpd00023_c0', 'cpd00033_c0', 'cpd00035_c0', 'cpd00039_c0', 'cpd00041_c0', 'cpd00051_c0', 'cpd00053_c0', 'cpd00054_c0', 'cpd00060_c0', 'cpd00065_c0', 'cpd00066_c0', 'cpd00069_c0', 'cpd00084_c0', 'cpd00107_c0', 'cpd00119_c0', 'cpd00129_c0', 'cpd00132_c0', 'cpd00156_c0', 'cpd00161_c0', 'cpd00322_c0', 'cpd00052_c0', 'cpd00038_c0', 'cpd00062_c0', 'cpd00115_c0', 'cpd00356_c0', 'cpd00241_c0', 'cpd00357_c0', 'cpd19001_c0', 'cpd30321_c0', 'cpd00163_c0', 'cpd16443_c0', 'cpd00080_c0', 'cpd01059_c0', 'cpd00604_c0', 'cpd00214_c0', 'cpd00536_c0', 'cpd01080_c0', 'cpd00104_c0', 'cpd00056_c0', 'cpd00016_c0', 'cpd00003_c0', 'cpd00004_c0', 'cpd00006_c0', 'cpd00005_c0', 'cpd00015_c0', 'cpd00050_c0', 'cpd00010_c0', 'cpd00087_c0', 'cpd02197_c0', 'cpd00201_c0', 'cpd00347_c0', 'cpd00125_c0', 'cpd00345_c0', 'cpd25914_c0', 'cpd16503_c0', 'cpd00017_c0', 'cpd00834_c0', 'cpd12844_d0', 'cpd00032_c0', 'cpd00130_c0', 'cpd00137_c0', 'cpd00331_c0', 'cpd00159_c0', 'cpd00205_c0', 'cpd00099_c0', 'cpd00012_c0', 'cpd00014_c0');
my %metabolites=(
"cpd00023_c0"=> -0.13,
"cpd00033_c0"=> -0.15,
"cpd00035_c0"=> -0.161,
"cpd00039_c0"=> -0.096,
"cpd00041_c0"=> -0.098,
"cpd00051_c0"=> -0.11,
"cpd00053_c0"=> -0.071,
"cpd00054_c0"=> -0.1,
"cpd00060_c0"=> -0.049,
"cpd00065_c0"=> -0.024,
"cpd00066_c0"=> -0.074,
"cpd00069_c0"=> -0.049,
"cpd00084_c0"=> -0.029,
"cpd00107_c0"=> -0.172,
"cpd00119_c0"=> -0.038,
"cpd00129_c0"=> -0.096,
"cpd00132_c0"=> -0.056,
"cpd00156_c0"=> -0.127,
"cpd00161_c0"=> -0.096,
"cpd00322_c0"=> -0.083,
"cpd00052_c0"=> -0.026,
"cpd00038_c0"=> -0.026,
"cpd00062_c0"=> -0.027,
"cpd00115_c0"=> -0.03,
"cpd00356_c0"=> -0.026,
"cpd00241_c0"=> -0.026,
"cpd00357_c0"=> -0.03,
"cpd19001_c0"=> -1.59,
"cpd30321_c0"=> -0.2271,
"cpd00163_c0"=> -0.766,
"cpd16443_c0"=> -0.0377,
"cpd00080_c0"=> -0.13,
"cpd01059_c0"=> -0.076,
"cpd00604_c0"=> -0.55,
"cpd00214_c0"=> -0.329,
"cpd00536_c0"=> -0.015,
"cpd01080_c0"=> -0.011,
"cpd00104_c0"=> -0.0001,
"cpd00056_c0"=> -0.0001,
"cpd00016_c0"=> -0.0001,
"cpd00003_c0"=> -0.0003,
"cpd00004_c0"=> -0.00015,
"cpd00006_c0"=> -0.00013,
"cpd00005_c0"=> -0.0001,
"cpd00015_c0"=> -0.0001,
"cpd00050_c0"=> -0.0001,
"cpd00010_c0"=> -0.000136,
"cpd00087_c0"=> -0.0001,
"cpd02197_c0"=> -0.0001,
"cpd00201_c0"=> -0.0001,
"cpd00347_c0"=> -0.0001,
"cpd00125_c0"=> -0.0001,
"cpd00345_c0"=> -0.0001,
"cpd25914_c0"=> -0.0001,
"cpd16503_c0"=> -0.000155,
"cpd00017_c0"=> -0.000543,
"cpd00834_c0"=> -0.001,
"cpd12844_d0"=> -0.000975,
"cpd00032_c0"=> -0.0757,
"cpd00130_c0"=> -0.037,
"cpd00137_c0"=> -0.013,
"cpd00331_c0"=> -0.086,
"cpd00159_c0"=> -0.039,
"cpd00205_c0"=> -0.307,
"cpd00099_c0"=> -0.21,
"cpd00012_c0"=> 0.218,
"cpd00014_c0"=> 0.766
);

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
	if ($cpd->{id} =~ m/_e\d+$/ || $cpd->{id} =~ m/cpd11416/ || $cpd->{id} =~ m/cpd02701/ ||  grep { $_ eq $cpd->{id}} keys %metabolites) {
		if ( grep { $_ eq $cpd->{id}} keys %metabolites) { #if biomass drain flux
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

#Adding flexible biomass drain constraints
my $drain_constraints;
foreach my $met (keys %metabolites){
	push(@{$problem->{constraints}},{
		name => "BCMax_".$met,
		type => "BiomassCompoundConstraint",
		variables => [$bio_obj->{bio1}->{variables}->{Flux},$cpd_obj->{$met}->{variables}->{Flux}],
		coefficients => [0.75 * abs($metabolites{$met}),-1],
		rhs => 0,
		sign => "<="
	});
	push(@{$problem->{constraints}},{
		name => "BCMin_".$met,
		type => "BiomassCompoundConstraint",
		variables => [$bio_obj->{bio1}->{variables}->{Flux},$cpd_obj->{$met}->{variables}->{Flux}],
		coefficients => [-0.75 * abs($metabolites{$met}),-1],
		rhs => 0,
		sign => ">="
	});
}

#Creating objective
my $rxnhash = {};
#Adding measured flux constraint
foreach my $rxnid (keys(%{$reaction_measured_hash})) {
	#Adding constraint: v - vm = vfit
	foreach my $dataset (keys(%{$datasets})) {
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
		foreach my $dataset (keys(%{$datasets})) {
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

#Adding kfit optimal flux constraint
foreach my $rxnid (keys(%{$reaction_protein_hash})) {
	if (!defined($rxnhash->{$rxnid})) {
		#Adding constraint: v - vm = vfit
		foreach my $dataset (keys(%{$datasets})) {
			my $full_id = $rxnid.$datasets->{$dataset};
			$rxnhash->{$full_id} = 1;
			#Adding kfit variable
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
			push(@{$problem->{objective}->{coefficients}},$vopt_coef);
			push(@{$problem->{objective}->{quadratic}},1);
		}
	}
}

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
$output->[6] .= "\n $drain_constraints";

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