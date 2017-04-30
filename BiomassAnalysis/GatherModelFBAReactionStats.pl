use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();

#my $rxndatafilename = "/Users/chenry/Dropbox/Exxon_Project_Data/Reactions.json";
#my $sortedpathways = "/Users/chenry/code/ResearchScripts/BiomassAnalysis/SortedPathways.txt";
#my $outputfile = "/Users/chenry/workspace/ModelRxnData.json";
#my $keggfile = "/Users/chenry/code/ModelSEEDDatabase/Pathways/KEGG.pathways";
my $rxndatafilename = "/homes/chenry/Reactions.json";
my $sortedpathways = "/homes/chenry/SortedPathways.txt";
my $outputfile = "/disks/p3dev3/ModelRxnData.json";
my $keggfile = "/homes/chenry/KEGG.pathways";
my $sortedpathways = Bio::KBase::ObjectAPI::utilities::LOADFILE($sortedpathways);
my $rxndata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($rxndatafilename)}));
my $full_rxn_data;
for (my $i=0; $i < @{$rxndata}; $i++) {
	$full_rxn_data->{$rxndata->[$i]->{id}} = $rxndata->[$i];
}
my $keggpathways = Bio::KBase::ObjectAPI::utilities::LOADFILE($keggfile);
my $pathwaytrans = {};
for (my $i=0; $i < @{$keggpathways}; $i++) {
	my $array = [split(/\t/,$keggpathways->[$i])];
	$array->[1] =~ s/map/rn/;
	$pathwaytrans->{$array->[1]} = $array->[2];
}

my $modelws = "chenry:narrative_1493181437626";
my $fbaws = "chenry:narrative_1493303324098";
my $media = "Carbon-D-Glucose";
my $mediaws = "KBaseMedia";

my $ws = Bio::KBase::kbaseenv::ws_client();
my $models = $ws->list_objects({
	workspaces => [$modelws],
	type => "KBaseFBA.FBAModel"
});

my $reaction_hash;
#for (my $i=0; $i < @{$models}; $i++) {
for (my $i=0; $i < 10; $i++) {
	print "Processing ".$models->[$i]->[1]."\n";
	#eval {
		my $model = $impl->util_get_object($modelws."/".$models->[$i]->[1]);
		my $rxns = $model->modelreactions();
		for (my $j=0; $j < @{$rxns}; $j++) {
			my $rxn = $rxns->[$j];
			my $id = $rxn->id();
			if ($id =~ m/(.+)_[a-z]\d+$/) {
				$id = $1;
			}
			if (!defined($reaction_hash->{$id})) {
				$reaction_hash->{$id} = {
					genomes => {},
					directions => {},
					equation => $rxn->equation(),
					definition => $rxn->definition(),
					pathway => undef,
					pathways => [],
					coupledbios => {},
					coupledrxns => {},
					mmclasses => {},
					comclasses => {},
					coupledrxn_pathways => {unassigned => 0},
					modelcount => 0,
					gapfillcount => 0,
					spontaneouscount => 0
				};
				if (defined($full_rxn_data->{$id})) {
					$reaction_hash->{$id}->{roles} = $full_rxn_data->{$id}->{roles}
				}
				if (defined($full_rxn_data->{$id}->{kegg_pathways})) {
					for (my $n=0; $n < @{$full_rxn_data->{$id}->{kegg_pathways}}; $n++) {
						push(@{$reaction_hash->{$id}->{pathways}},$pathwaytrans->{$full_rxn_data->{$id}->{kegg_pathways}->[$n]});
					}
					for (my $m=0; $m < @{$sortedpathways}; $m++) {
						for (my $n=0; $n < @{$full_rxn_data->{$id}->{kegg_pathways}}; $n++) {
							if ($pathwaytrans->{$full_rxn_data->{$id}->{kegg_pathways}->[$n]} eq $sortedpathways->[$m]) {
								$reaction_hash->{$id}->{pathway} = $sortedpathways->[$m];
								last;
							}
						}
						if (defined($reaction_hash->{$id}->{pathway})) {
							last;
						}
					}
				}
			}
			if (!defined($reaction_hash->{$id}->{directions}->{$rxn->direction()})) {
				$reaction_hash->{$id}->{directions}->{$rxn->direction()} = 0;
			}
			$reaction_hash->{$id}->{directions}->{$rxn->direction()}++;
			$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]} = {
				direction => $rxn->direction(),
				gpr => $rxn->gprString(),
				genes => $rxn->featureIDs(),
				gapfill => $rxn->gapfillString(),
				mmclass => undef,
				comclass => undef,
				mmflux => undef,
				comflux => undef,
				biomassdep => [],
				coupledrxn => []
			};
			$reaction_hash->{$id}->{modelcount}++;
			if ($rxn->gprString() eq "Unknown") {
				if ($rxn->gapfillString() eq "") {
					$reaction_hash->{$id}->{spontaneouscount}++;
				} else {
					$reaction_hash->{$id}->{gapfillcount}++;
				}
			}
		}
		my $suffixes = ["sensfba","mmfva","comfva"];
		for (my $j=0; $j < @{$suffixes}; $j++) {
			my $fba = $impl->util_get_object($fbaws."/".$models->[$i]->[1].".".$suffixes->[$j]);
			$rxns = $fba->FBAReactionVariables();
			for (my $k=0; $k < @{$rxns}; $k++) {
				my $id = $rxns->[$k]->modelreaction()->id();
				if ($id =~ m/(.+)_[a-z]\d+$/) {
					$id = $1;
				}
				if ($j == 0) {
					if (defined($rxns->[$k]->biomass_dependencies())) {
						for (my $m=0; $m < @{$rxns->[$k]->biomass_dependencies()}; $m++) {
							push(@{$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{biomassdep}},$rxns->[$k]->biomass_dependencies()->[$m]->[1]);
							if (!defined($reaction_hash->{$id}->{coupledbios}->{$rxns->[$k]->biomass_dependencies()->[$m]->[1]})) {
								$reaction_hash->{$id}->{coupledbios}->{$rxns->[$k]->biomass_dependencies()->[$m]->[1]} = [0,$model->getObject("modelcompounds",$rxns->[$k]->biomass_dependencies()->[$m]->[1])->name()];
							}
							$reaction_hash->{$id}->{coupledbios}->{$rxns->[$k]->biomass_dependencies()->[$m]->[1]}->[0]++;
						}
					}
					if (defined($rxns->[$k]->coupled_reactions())) {
						for (my $m=0; $m < @{$rxns->[$k]->coupled_reactions()}; $m++) {
							my $cid = $rxns->[$k]->coupled_reactions()->[$m];
							if ($cid =~ m/(.+)_[a-z]\d+$/) {
								$cid = $1;
							}	
							push(@{$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{coupledrxn}},$cid);
							if (!defined($reaction_hash->{$id}->{coupledrxns}->{$cid})) {
								$reaction_hash->{$id}->{coupledrxns}->{$cid} = 0;
							}
							$reaction_hash->{$id}->{coupledrxns}->{$cid}++;
							if (defined($reaction_hash->{$cid}->{pathway})) {
								if (!defined($reaction_hash->{$id}->{coupledrxn_pathways}->{$reaction_hash->{$cid}->{pathway}})) {
									$reaction_hash->{$id}->{coupledrxn_pathways}->{$reaction_hash->{$cid}->{pathway}} = 0;
								}
								$reaction_hash->{$id}->{coupledrxn_pathways}->{$reaction_hash->{$cid}->{pathway}}++;
							} else {
								$reaction_hash->{$id}->{coupledrxn_pathways}->{unassigned}++;
							}
						}
					}
				} elsif ($j == 1) {
					$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{mmclass} = $rxns->[$k]->class();
					$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{mmflux} = $rxns->[$k]->value();
					if (!defined($reaction_hash->{$id}->{mmclasses}->{$rxns->[$k]->class()})) {
						$reaction_hash->{$id}->{mmclasses}->{$rxns->[$k]->class()} = 0;
					}
					$reaction_hash->{$id}->{mmclasses}->{$rxns->[$k]->class()}++;
				} elsif ($j == 2) {
					$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{comclass} = $rxns->[$k]->class();
					$reaction_hash->{$id}->{genomes}->{$models->[$i]->[1]}->{comflux} = $rxns->[$k]->value();
					if (!defined($reaction_hash->{$id}->{comclasses}->{$rxns->[$k]->class()})) {
						$reaction_hash->{$id}->{comclasses}->{$rxns->[$k]->class()} = 0;
					}
					$reaction_hash->{$id}->{comclasses}->{$rxns->[$k]->class()}++;
				}
			}
		}
	#};
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE($outputfile,[Bio::KBase::ObjectAPI::utilities::TOJSON($reaction_hash,1)]);