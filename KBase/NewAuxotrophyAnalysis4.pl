#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $params = {
	workspace => "weise:narrative_1531417693962",
	genome_workspace => "weise:narrative_1531417693962",
	genome_ids => []
};

my $output = Bio::KBase::kbaseenv::list_objects({
	workspaces => ["weise:narrative_1531417693962"]
});

for (my $i=0; $i < @{$output}; $i++) {
	if ($output->[$i]->[2] =~ m/KBaseGenomes\.Genome/) {
		push(@{$params->{genome_ids}},$output->[$i]->[1]);
	}
}

my $corerxns = [
["rxn05319","="],
["rxn05467","="],
["rxn05468",">"],
["rxn00011","<"],
["rxn00077","="],
["rxn00083","="],
["rxn00102","="],
["rxn00122",">"],
["rxn00147",">"],
["rxn00148","<"],
["rxn00151",">"],
["rxn00157","<"],
["rxn00159",">"],
["rxn00160",">"],
["rxn00161",">"],
["rxn00162",">"],
["rxn00171","="],
["rxn00172","<"],
["rxn00173",">"],
["rxn00175",">"],
["rxn00178",">"],
["rxn00199",">"],
["rxn00216","="],
["rxn00224","="],
["rxn00225","<"],
["rxn00247",">"],
["rxn00248","="],
["rxn00250","="],
["rxn00251","="],
["rxn00256","<"],
["rxn00265",">"],
["rxn00285","="],
["rxn00288",">"],
["rxn00305",">"],
["rxn00324","<"],
["rxn00330","<"],
["rxn00336",">"],
["rxn00371",">"],
["rxn00392","="],
["rxn00441",">"],
["rxn00459","="],
["rxn00499","="],
["rxn00500","="],
["rxn00505",">"],
["rxn00506",">"],
["rxn00512","<"],
["rxn00543","="],
["rxn00544","<"],
["rxn00545",">"],
["rxn00548","="],
["rxn00549",">"],
["rxn00558","="],
["rxn00568","<"],
["rxn00569","<"],
["rxn00602","<"],
["rxn00604",">"],
["rxn00615",">"],
["rxn00616","="],
["rxn00623","="],
["rxn00670","="],
["rxn00747","="],
["rxn00762","="],
["rxn00763","="],
["rxn00770","="],
["rxn00777","="],
["rxn00779",">"],
["rxn00781","="],
["rxn00782","<"],
["rxn00785","="],
["rxn00786","="],
["rxn00799","="],
["rxn00868","<"],
["rxn00869","<"],
["rxn00871",">"],
["rxn00874","="],
["rxn00875",">"],
["rxn00935","="],
["rxn00973",">"],
["rxn00974","="],
["rxn00985","="],
["rxn00990",">"],
["rxn00994","<"],
["rxn00995",">"],
["rxn01100","="],
["rxn01106","="],
["rxn01115",">"],
["rxn01116","="],
["rxn01121",">"],
["rxn01123",">"],
["rxn01130",">"],
["rxn01169","="],
["rxn01187",">"],
["rxn01200","="],
["rxn01236","<"],
["rxn01241","="],
["rxn01275",">"],
["rxn01333","="],
["rxn01334","="],
["rxn01343",">"],
["rxn01387","="],
["rxn01388","="],
["rxn01452","<"],
["rxn01476",">"],
["rxn01477",">"],
["rxn01480","="],
["rxn01806",">"],
["rxn01870","="],
["rxn01871","<"],
["rxn01872","="],
["rxn01975",">"],
["rxn01977","="],
["rxn01996","="],
["rxn02112","="],
["rxn02113",">"],
["rxn02167",">"],
["rxn02168","="],
["rxn02185","="],
["rxn02314","="],
["rxn02342",">"],
["rxn02359","="],
["rxn02376","="],
["rxn02380","="],
["rxn02527",">"],
["rxn03240","="],
["rxn03249","="],
["rxn03250","="],
["rxn03643","="],
["rxn03644","="],
["rxn03884",">"],
["rxn03978","="],
["rxn04713","="],
["rxn04794","="],
["rxn05040",">"],
["rxn05145",">"],
["rxn05209","="],
["rxn05226",">"],
["rxn05313",">"],
["rxn05528",">"],
["rxn05559","<"],
["rxn05561",">"],
["rxn05573",">"],
["rxn05581","="],
["rxn05602","="],
["rxn05625","="],
["rxn05627",">"],
["rxn05937","="],
["rxn05938","<"],
["rxn05939","="],
["rxn06493","="],
["rxn06526","="],
["rxn06777","="],
["rxn07191","="],
["rxn08094",">"],
["rxn08178",">"],
["rxn08179",">"],
["rxn08527","="],
["rxn08556",">"],
["rxn08557",">"],
["rxn08655",">"],
["rxn08656",">"],
["rxn08901",">"],
["rxn08971",">"],
["rxn08975",">"],
["rxn08977",">"],
["rxn09272",">"],
["rxn10042","="],
["rxn10114",">"],
["rxn10118",">"],
["rxn10120",">"],
["rxn10121",">"],
["rxn10122",">"],
["rxn10125",">"],
["rxn10126",">"],
["rxn10127",">"],
["rxn10154",">"],
["rxn11937","<"],
#["rxn14412",">"],
#["rxn14414",">"],
#["rxn14416",">"],
#["rxn14418",">"],
#["rxn14419",">"],
#["rxn14420",">"],
#["rxn14421",">"],
#["rxn14422",">"],
#["rxn14423",">"],
#["rxn14424",">"],
#["rxn14425",">"],
#["rxn14426",">"],
#["rxn14427",">"],
#["rxn14428",">"]
];

my $biocpds = [
{modelcompound_ref => "~/modelcompounds/id/cpd00022_c0",coefficient => -3.7478},
{modelcompound_ref => "~/modelcompounds/id/cpd00003_c0",coefficient => -3.547},
{modelcompound_ref => "~/modelcompounds/id/cpd00236_c0",coefficient => -0.8977},
{modelcompound_ref => "~/modelcompounds/id/cpd00010_c0",coefficient => 3.7478},
{modelcompound_ref => "~/modelcompounds/id/cpd00101_c0",coefficient => -0.8977},
{modelcompound_ref => "~/modelcompounds/id/cpd00061_c0",coefficient => -0.5191},
{modelcompound_ref => "~/modelcompounds/id/cpd00079_c0",coefficient => -0.205},
{modelcompound_ref => "~/modelcompounds/id/cpd00032_c0",coefficient => -1.7867},
{modelcompound_ref => "~/modelcompounds/id/cpd00004_c0",coefficient => 3.547},
{modelcompound_ref => "~/modelcompounds/id/cpd00072_c0",coefficient => -0.0709},
{modelcompound_ref => "~/modelcompounds/id/cpd00005_c0",coefficient => -1.8225},
{modelcompound_ref => "~/modelcompounds/id/cpd00002_c0",coefficient => -41.257},
{modelcompound_ref => "~/modelcompounds/id/cpd00009_c0",coefficient => 41.257},
{modelcompound_ref => "~/modelcompounds/id/cpd00008_c0",coefficient => 41.257},
{modelcompound_ref => "~/modelcompounds/id/cpd11416_c0",coefficient => 1},
{modelcompound_ref => "~/modelcompounds/id/cpd00067_c0",coefficient => -35.8875},
{modelcompound_ref => "~/modelcompounds/id/cpd00102_c0",coefficient => -0.129},
{modelcompound_ref => "~/modelcompounds/id/cpd00169_c0",coefficient => -1.496},
{modelcompound_ref => "~/modelcompounds/id/cpd00006_c0",coefficient => 1.8225},
{modelcompound_ref => "~/modelcompounds/id/cpd00024_c0",coefficient => -1.0789},
{modelcompound_ref => "~/modelcompounds/id/cpd00020_c0",coefficient => -2.8328},
{modelcompound_ref => "~/modelcompounds/id/cpd00001_c0",coefficient => 41.257}
];

my $media = $handler->util_get_object(Bio::KBase::utilities::conf("ModelSEED","default_media_workspace")."/Carbon-D-Glucose");
my $genomes = $params->{genome_ids};
my $transporthash = {};
my $cpddatahash = Bio::KBase::utilities::compound_hash();
my $auxotrophy_threshold_hash = Bio::KBase::constants::auxotrophy_thresholds();
my $cofarray = Bio::KBase::constants::cofactors();
my $cofactors = {};
foreach my $id (@{$cofarray}) {
	$cofactors->{$id."_c0"} = 1;
}

my $auxotrophy_hash;
my $template_trans = Bio::KBase::constants::template_trans();
for (my $i=0; $i < @{$genomes}; $i++) {
	my $datachannel = {};
	my $genomeobj = $handler->util_get_object(Bio::KBase::utilities::buildref($genomes->[$i],$params->{genome_workspace}));
	my $genomeid = $genomeobj->_wsname();
	$genomes->[$i] = $genomeid;
	print "Processing ".$genomeid."\n";
	my $current_media = $media->cloneObject();
	my $tid = $template_trans->{$genomeobj->template_classification()};
	#Building draft model
	Bio::KBase::ObjectAPI::functions::func_build_metabolic_model({
		template_id => $tid,
		template_workspace => Bio::KBase::utilities::conf("ModelSEED","default_template_workspace"),
		workspace => "NULL",
		fbamodel_output_id => $genomeid.".fbamodel",
		genome_id => $genomeid,
		genome_workspace => $params->{genome_workspace},
		media_id => "Carbon-D-Glucose",
		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
		gapfill_model => 0
	},$datachannel);
	#Adding biotin to model compounds
	if (!defined($datachannel->{fbamodel}->queryObject("modelcompounds",{id => "cpd00104_c0"}))) {
		my $cpddata = {
			id => "cpd00104_c0",
			modelcompartment_ref => "~/modelcompartments/id/c0",
			compound_ref => "~/template_ref/compcompounds/id/cpd00104_c"
		};
		my $list = ["charge","formula","inchikey","smiles","name"];
		for (my $j=0; $j < @{$list}; $j++) {
			if (defined($cpddatahash->{cpd00104}->{$list->[$j]})) {
				$cpddata->{$list->[$j]} = $cpddatahash->{cpd00104}->{$list->[$j]};
			}
		}
		$datachannel->{fbamodel}->add("modelcompounds",$cpddata);
	}
	#Adding biotin to primary biomass
	$datachannel->{fbamodel}->biomasses()->[0]->add("biomasscompounds",{
		modelcompound_ref => "~/modelcompounds/id/cpd00104_c0",
		coefficient => -0.00001
	});
#	#Adding core biomass
#	$datachannel->{fbamodel}->add("biomasses",{
#		id => "bio2",
#		name => "bio2",
#		other => 1,
#		dna => 0,
#		rna => 0,
#		protein => 0,
#		cellwall => 0,
#		lipid => 0,
#		cofactor => 0,
#		energy => 0,
#		biomasscompounds => $biocpds,
#		removedcompounds => []
#	});
#	#Adding core reactions
#	foreach my $rxn (@{$corerxns}) {
#		if (!defined($datachannel->{fbamodel}->queryObject("modelreactions",{id => $rxn->[0]."_c0"}))) {
#			$datachannel->{fbamodel}->addModelReaction({
#				reaction => $rxn->[0],
#				direction => $rxn->[1],
#				addReaction => 1
#			});
#		}
#	}
#	#Gapfill minimal media with core biomass
#	Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({
#		target_reaction => "bio2",
#		minimum_target_flux => 0.1,
#		number_of_solutions => 1,
#		workspace => "NULL",
#		fbamodel_id => $genomeid.".fbamodel",
#		fbamodel_output_id => $genomeid.".fbamodel",
#		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
#		media_id => "Carbon-D-Glucose"
#	},$datachannel->{fbamodel});
	#Gapfill minimal media with genome-scale biomass
	Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({
		target_reaction => "bio1",
		minimum_target_flux => 0.1,
		number_of_solutions => 1,
		workspace => "NULL",
		fbamodel_id => $genomeid.".fbamodel",
		fbamodel_output_id => $genomeid.".fbamodel",
		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
		media_id => "Carbon-D-Glucose"
	},$datachannel->{fbamodel});
	#Adding transport reactions
	foreach my $rxn (keys(%{$transporthash})) {
		if (!defined($datachannel->{fbamodel}->queryObject("modelreactions",{id => $rxn."_c0"}))) {
			$datachannel->{fbamodel}->addModelReaction({
				reaction => $rxn,
				direction => "=",
				addReaction => 1
			});
		}
	}
	#Populating reaction and compound hashes
	my $rxns = $datachannel->{fbamodel}->modelreactions();
	my $rxnhash;
	for (my $j=0; $j < @{$rxns}; $j++) {
		$rxnhash->{$rxns->[$j]->id()} = $rxns->[$j];
	}
	my $cpds = $datachannel->{fbamodel}->modelcompounds();
	my $cpdhash;
	for (my $j=0; $j < @{$cpds}; $j++) {
		$cpdhash->{$cpds->[$j]->id()} = $cpds->[$j];
	}
	#Running auxotrophy prediction flux balance analysis
	my $filelist = [];
	Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
		workspace => "NULL",
		fbamodel_id => $genomeid.".fbamodel",
		fba_output_id => $genomeid.".fba_auxo1",
		media_id => "Carbon-D-Glucose",
		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
		target_reaction => "bio1",
		metabolite_production_analysis => 1,
		source_metabolite_list => ["cpd00103_c0","cpd00171_c0","cpd00146_c0","cpd00020_c0","cpd00024_c0","cpd00169_c0","cpd00102_c0","cpd00072_c0","cpd00032_c0",
			"cpd00079_c0","cpd00022_c0","cpd00236_c0","cpd00101_c0","cpd00061_c0","cpd00041_c0","cpd00002_c0","cpd00038_c0","cpd00023_c0","cpd00053_c0"],
		target_metabolite_list => ["cpd00017_c0","cpd00033_c0","cpd00054_c0","cpd00161_c0","cpd00084_c0","cpd00119_c0","cpd00060_c0","cpd00051_c0","cpd00129_c0",
			"cpd00118_c0","cpd00132_c0","cpd00016_c0","cpd00056_c0","cpd00220_c0","cpd00028_c0","cpd00166_c0","cpd00557_c0","cpd00039_c0","cpd00069_c0",
			"cpd00066_c0","cpd00065_c0","cpd00393_c0","cpd00156_c0","cpd00107_c0"],
	},$datachannel);
	$filelist->[0] = $datachannel->{fba}->outputfiles()->{MetaboliteProductionResults};
	Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
		workspace => "NULL",
		fbamodel_id => $genomeid.".fbamodel",
		fba_output_id => $genomeid.".fba_auxo1",
		media_id => "Carbon-D-Glucose",
		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
		target_reaction => "bio1",
		metabolite_production_analysis => 1,
		source_metabolite_list => ["cpd00103_c0","cpd00171_c0","cpd00146_c0","cpd00020_c0","cpd00024_c0","cpd00169_c0","cpd00102_c0","cpd00072_c0","cpd00032_c0",
			"cpd00079_c0","cpd00022_c0","cpd00236_c0","cpd00101_c0","cpd00061_c0","cpd00041_c0","cpd00002_c0","cpd00038_c0","cpd00023_c0","cpd00053_c0","cpd00054_c0"],
		target_metabolite_list => ["cpd00065"],
	},$datachannel);
	$filelist->[1] = $datachannel->{fba}->outputfiles()->{MetaboliteProductionResults};
	Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
		workspace => "NULL",
		fbamodel_id => $genomeid.".fbamodel",
		fba_output_id => $genomeid.".fba_auxo2",
		media_id => "Carbon-D-Glucose",
		media_workspace => Bio::KBase::utilities::conf("ModelSEED","default_media_workspace"),
		target_reaction => "bio1",
		metabolite_production_analysis => 1,
		source_metabolite_list => ["cpd00800_c0","cpd00103_c0","cpd00171_c0","cpd00146_c0","cpd00020_c0","cpd00024_c0","cpd00169_c0","cpd00102_c0","cpd00072_c0","cpd00032_c0",
		"cpd00079_c0","cpd00022_c0","cpd00236_c0","cpd00101_c0","cpd00061_c0","cpd00041_c0","cpd00002_c0","cpd00038_c0","cpd00023_c0","cpd00053_c0",
		"cpd00017_c0","cpd00051_c0","cpd00084_c0","cpd00065_c0","cpd00161_c0","cpd00156_c0","cpd00800_c0","cpd00054_c0"],
		target_metabolite_list => ["cpd00065","cpd00644_c0","cpd00264_c0","cpd00042_c0","cpd00003_c0","cpd00104_c0","cpd00322_c0"],
	},$datachannel);
	$filelist->[2] = $datachannel->{fba}->outputfiles()->{MetaboliteProductionResults};
	my $auxotrophy_results = {};
	my $translation = {
		cpd00016 => "cpd00215"
	};
	for (my $j=0; $j < @{$filelist}; $j++) {
		for (my $k=1; $k < @{$filelist->[$j]}; $k++) {
			my $array = [split(/\t/,$filelist->[$j]->[$k])];
			my $cpd = $array->[1];
			$cpd =~ s/_c0//;
			if (defined($translation->{$cpd})) {
				$cpd = $translation->{$cpd};
			}
			if ($array->[2] ne "none") {
				$auxotrophy_results->{$cpd} = {
					totalrxn => 0,
					gfrxn => 0,
					rxns => {},
					auxotrophic => 0,
					cofactor_reactants => [],
					main_reactants => [],
					cofactor_products => [],
					main_products => [],
					name => $cpddatahash->{$cpd}->{name}
				};
				my $rxns = [split(/;/,$array->[2])];
				my $rxncpdhash = {};
				for (my $m=0; $m < @{$rxns}; $m++) {
					if ($rxns->[$m] =~ m/(.)(.+)(_[a-zA-Z]\d):(.+)/) {
						my $rxnid = $2;
						my $comp = $3;
						$auxotrophy_results->{$cpd}->{rxns}->{$rxnid} = {
							flux => $4,
							multiplier => 1,
							cofactor_reactants => [],
							main_reactants => [],
							cofactor_products => [],
							main_products => [],
							gfrxn => 0
						};
						if ($auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{flux} < 0) {
							$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier} = -1;
						}
						$auxotrophy_results->{$cpd}->{totalrxn}++;
						if ($rxnhash->{$rxnid.$comp}->gprString() eq "Unknown") {
							$auxotrophy_results->{$cpd}->{gfrxn}++;
							$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{gfrxn} = 1;
						}
						my $rgts = $rxnhash->{$rxnid.$comp}->modelReactionReagents();
						#Populating data for reaction stoichiometry
						for (my $n=0; $n < @{$rgts}; $n++) {
							if (!defined($rxncpdhash->{$rgts->[$n]->modelcompound()->id()})) {
								$rxncpdhash->{$rgts->[$n]->modelcompound()->id()} = 0;
							}
							$rxncpdhash->{$rgts->[$n]->modelcompound()->id()} += $rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{flux};
							if (defined($cofactors->{$cpd})) {
								if ($rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier} < 0) {
									push(@{$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{cofactor_reactants}},[$rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier},$cpd,$cpddatahash->{$cpd}->{name}]);
								} else {
									push(@{$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{cofactor_products}},[$rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier},$cpd,$cpddatahash->{$cpd}->{name}]);
								}
							} else {
								if ($rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier} < 0) {
									push(@{$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{main_reactants}},[$rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier},$cpd,$cpddatahash->{$cpd}->{name}]);
								} else {
									push(@{$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{main_products}},[$rgts->[$n]->coefficient()*$auxotrophy_results->{$cpd}->{rxns}->{$rxnid}->{multiplier},$cpd,$cpddatahash->{$cpd}->{name}]);
								}
							}
						}
					}
				}
				#Populating data for overall reaction for metabolite
				foreach my $cpd (keys(%{$rxncpdhash})) {
					if (defined($cofactors->{$cpd})) {
						if ($rxncpdhash->{$cpd} < 0) {
							push(@{$auxotrophy_results->{$cpd}->{cofactor_reactants}},[$rxncpdhash->{$cpd},$cpd,$cpddatahash->{$cpd}->{name}]);
						} else {
							push(@{$auxotrophy_results->{$cpd}->{cofactor_products}},[$rxncpdhash->{$cpd},$cpd,$cpddatahash->{$cpd}->{name}]);
						}
					} else {
						if ($rxncpdhash->{$cpd} < 0) {
							push(@{$auxotrophy_results->{$cpd}->{main_reactants}},[$rxncpdhash->{$cpd},$cpd,$cpddatahash->{$cpd}->{name}]);
						} else {
							push(@{$auxotrophy_results->{$cpd}->{main_products}},[$rxncpdhash->{$cpd},$cpd,$cpddatahash->{$cpd}->{name}]);
						}
					}
				}
				if ($auxotrophy_results->{$cpd}->{gfrxn} >= $auxotrophy_threshold_hash->{$cpd}->[1]) {
					$auxotrophy_results->{$cpd}->{auxotrophic} = 1;
					$current_media->add("mediacompounds",{
						compound_ref => "kbase/default/compounds/id/".$cpd,
						id => $cpd,
						name => $cpddatahash->{$cpd}->{name},
						concentration => 0.001,
						maxFlux => 100,
						minFlux => -100
					});
				}
			}
			
		}
	}
	$current_media->parent($handler->util_store());
	my $wsmeta = $handler->util_save_object($datachannel->{fbamodel},$params->{workspace}."/".$genomeid.".fbamodel");
	$datachannel->{fba}->{fbamodel_ref} = $params->{workspace}."/".$genomeid.".fbamodel";
	$wsmeta = $handler->util_save_object($datachannel->{fba},$params->{workspace}."/".$genomeid.".fba");
	$wsmeta = $handler->util_save_object($current_media,$params->{workspace}."/".$genomeid.".auxo_media");
	$auxotrophy_hash->{$genomeid} = $auxotrophy_results;
}

#print TOJSONBio::KBase::ObjectAPI::utilities::TOJSON($auxotrophy_hash,1);

my $cpdlist = [keys(%{$auxotrophy_threshold_hash})];
print "Compound\tName\t".join("\t",@{$genomes})."\n";
for (my $i=0; $i < @{$cpdlist}; $i++) {
	print $cpdlist->[$i]."\t".$cpddatahash->{$cpdlist->[$i]}->{name};
	for (my $j=0; $j < @{$genomes}; $j++) {
		if (defined($auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]})) {
			my $rxn = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{totalrxn};
			my $gf = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{gfrxn};
			my $auxo = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{auxotrophic};
			print "\t".$rxn;
		} else {
			print "\t--";
		}
	}
	print "\n";
}
for (my $i=0; $i < @{$cpdlist}; $i++) {
	print $cpdlist->[$i]."\t".$cpddatahash->{$cpdlist->[$i]}->{name};
	for (my $j=0; $j < @{$genomes}; $j++) {
		if (defined($auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]})) {
			my $rxn = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{totalrxn};
			my $gf = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{gfrxn};
			my $auxo = $auxotrophy_hash->{$genomes->[$j]}->{$cpdlist->[$i]}->{auxotrophic};
			print "\t".$gf;
		} else {
			print "\t--";
		}
	}
	print "\n";
}

#Auxotrophy analysis with core biomass as objective
#print "Class\tCompound\tAve gf\n";
#for (my $i=0; $i < @{$genomes}; $i++) {
#	print "\t".$genomes->[$i];
#}
#print "\n";
#for (my $i=0; $i < @{$cpddata}; $i++) {
#	print $cpddata->[$i]->{class}."\t".$cpddata->[$i]->{name}." (".$cpddata->[$i]->{id}.")\t".$cpddata->[$i]->{avegf};
#	for (my $j=0; $j < @{$genomes}; $j++) {
#		if (defined($auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]})) {
#			print "\t".$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{gfrxn}."/".$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{rxn}."/".$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{auxotrophic};
#		} else {
#			print "\t-";
#		}
#	}
#	print "\n";
#}
#
#my $htmlreport = "<html><head><script type='text/javascript' src='https://www.google.com/jsapi'></script><script type='text/javascript'>google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});google.setOnLoadCallback(drawDashboard);";
#$htmlreport .= "function drawDashboard() {var data = new google.visualization.DataTable();";
#$htmlreport .= "data.addColumn('string','Class');";
#$htmlreport .= "data.addColumn('string','Essential metabolite');";
#for (my $i=0; $i < @{$genomes}; $i++) {
#	$htmlreport .= "data.addColumn('number','".$genomes->[$i]."');";
#}
#$htmlreport .= "data.addRows([";
#for (my $i=0; $i < @{$cpddata}; $i++) {
#	my $row = [];
#	$htmlreport .= '["'.$cpddata->[$i]->{class}.'","'.$cpddata->[$i]->{name}." (".$cpddata->[$i]->{id}.")".'",';
#	for (my $j=0; $j < @{$genomes}; $j++) {
#		if (defined($auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]})) {
#			if ($auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{auxotrophic} == 1) {
#				push(@{$row},'{v:1,f:"Gapfilling: '.$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{gfrxn}.'<br>Reactions: '.$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{rxn}.'<br>Auxotrophic"}');	
#			} else {
#				push(@{$row},'{v:0,f:"Gapfilling: '.$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{gfrxn}.'<br>Reactions: '.$auxotrophy_hash->{$cpddata->[$i]->{id}}->{$genomes->[$j]}->{rxn}.'"}');
#			}
#		} else {
#			push(@{$row},'{v:0,f:""}');
#		}
#	}
#	$htmlreport .= join(',',@{$row}).'],';
#}
#$htmlreport .= "]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});";
#$htmlreport .= "var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});";
#$htmlreport .= "var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});";
#$htmlreport .= "var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');";
#for (my $j=0; $j < @{$genomes}; $j++) {
#	$htmlreport .= "formatter.format(data, ".($j+2).");";
#}
#$htmlreport .= "dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>";
#$htmlreport .= "<body><h4>Results from auxotrophy prediction on all genomes</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>";
#print $htmlreport;
