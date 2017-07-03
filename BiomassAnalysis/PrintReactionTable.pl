use strict;
use fba_tools::fba_toolsImpl;
use Data::Dumper;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();

my $directory = "/Users/janakaanl/MSRepo/ResearchScripts/IGap/";
my $cpddata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($directory."Compounds.json")}));
my $globalcpdhash;
for (my $i=0; $i < @{$cpddata}; $i++) {
	$globalcpdhash->{$cpddata->[$i]->{id}} = $cpddata->[$i];
	#print &Dumper ($globalcpdhash->{$cpddata->[$i]->{id}})."\n";
}
my $rxnmodeldata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($directory."ModelFBAData.json")}));


my $pathway_file = Bio::KBase::ObjectAPI::utilities::LOADFILE($directory."RankedKEGGPathways.txt");

my $pathwayhash = {};
my $pathwayarray = [];
my $pathwayid = {};
my $forced_pathway_assignments = {};
for (my $i=1; $i < @{$pathway_file}; $i++) {
	my $array = [split(/\t/,$pathway_file->[$i])];
	push(@{$pathwayarray},$array->[1]);
	$pathwayhash->{$array->[0]} = $array->[1];
	$pathwayid->{$array->[1]} = $array->[0];
	for (my $j=3; $j < @{$array}; $j++) {
		$pathwayhash->{$array->[$j]} = $array->[1];
	}
}
my $rxndata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($directory."Reactions.json")}));


my $rxndatahash = {};
for (my $i=0; $i < @{$rxndata}; $i++) {
	$rxndatahash->{$rxndata->[$i]->{id}} = $rxndata->[$i];
}
#Print the reactions, organized and prioritized by pathway
open(my $fout,">",$directory."ModelReactions.txt");
print $fout "Pathway\tReaction\tEquation\tReversibility\tStatus\tKEGG ID\tMetaCyc ID\tName\tEC numbers\tRoles";
print  "Pathway\tReaction\tEquation\tReversibility\tStatus\tKEGG ID\tMetaCyc ID\tName\tEC numbers\tRoles";

my $initial_genomes = Bio::KBase::kbaseenv::ws_client()->list_objects({
       	workspaces => ["janakakbase:narrative_1499012973666"],
        type => "KBaseGenomes.Genome",
});

my $genomelist;
my $genome_hash = {};
for (my $i=0; $i < @{$initial_genomes}; $i++) {
        if ($initial_genomes->[$i]->[1] =~ m/(.+)\.RAST$/) {
                $genome_hash->{$initial_genomes->[$i]->[1]} = $1;
                push(@{$genomelist},$initial_genomes->[$i]->[1]);
        }
}

my $pathwaymodelhash;
for (my $i=0; $i < @{$genomelist}; $i++) {
	print $fout "\t".$genomelist->[$i];
	#print "\t".$genomelist->[$i];
}
print $fout "\n";
#print "\n";
#die;
my $donerxn = {};
for (my $i=0; $i < @{$pathwayarray}; $i++) {
	 my $list = [];
	 foreach my $rxn (keys(%{$rxnmodeldata})) {
	 	$rxn =~ s/_c0//;
	 	if (!defined($donerxn->{$rxn})) {
	 		if (defined($forced_pathway_assignments->{$rxn})) {
	 			if ($forced_pathway_assignments->{$rxn} eq $pathwayarray->[$i]) {
	 				$donerxn->{$rxn} = 1;
					(my $line,my $priority) = &print_reaction($rxn,$pathwayarray->[$i]);
					push(@{$list},{line => $line,priority => $priority});
	 			}
	 		} elsif (defined($rxndatahash->{$rxn})) {
				if (defined($rxndatahash->{$rxn}->{kegg_pathways})) {
					for (my $j=0; $j < @{$rxndatahash->{$rxn}->{kegg_pathways}}; $j++) {
						if (defined($pathwayhash->{$rxndatahash->{$rxn}->{kegg_pathways}->[$j]})) {
							if ($pathwayhash->{$rxndatahash->{$rxn}->{kegg_pathways}->[$j]} eq $pathwayarray->[$i]) {
								$donerxn->{$rxn} = 1;
								(my $line,my $priority) = &print_reaction($rxn,$pathwayarray->[$i]);
								push(@{$list},{line => $line,priority => $priority});;
							}
						}
					}
				} elsif (defined($rxndatahash->{$rxn}->{metacyc_pathways})) {
					for (my $j=0; $j < @{$rxndatahash->{$rxn}->{metacyc_pathways}}; $j++) {
						if (defined($pathwayhash->{$rxndatahash->{$rxn}->{metacyc_pathways}->[$j]})) {
							if ($pathwayhash->{$rxndatahash->{$rxn}->{metacyc_pathways}->[$j]} eq $pathwayarray->[$i]) {
								$donerxn->{$rxn} = 1;
								(my $line,my $priority) = &print_reaction($rxn,$pathwayarray->[$i]);
								push(@{$list},{line => $line,priority => $priority});
							}
						}
					}
				}
			}
	 	}
	}
	my $sortedlines = [sort { $b->{priority} <=> $a->{priority} } @{$list} ];
	for (my $j=0; $j < @{$sortedlines}; $j++) {
		print $fout $sortedlines->[$j]->{line}."\n";
	}
}
my $list;
foreach my $rxn (keys(%{$rxnmodeldata})) {
	$rxn =~ s/_c0//;
	if (!defined($donerxn->{$rxn})) {
		(my $line,my $priority) = &print_reaction($rxn);
		push(@{$list},{line => $line,priority => $priority});
	}
}
my $sortedlines = [sort { $b->{priority} <=> $a->{priority} } @{$list} ];
for (my $i=0; $i < @{$sortedlines}; $i++) {
	print $fout $sortedlines->[$i]->{line}."\n";
}
close($fout);

sub print_reaction {
	my ($id,$pathway) = @_;
	#Printing general reaction data
	my $line = "";
	if (defined($pathway)) {
		$line .= $pathway."\t".$id."\t";
	} else {
		$line .= "Unassigned\t".$id."\t";
	}
	if (defined($rxndatahash->{$id})) {
		my $rxnobj = $rxndatahash->{$id};
		my $headings = [
			["definition",0],
			["reversibility",0],
			["status",0],
			["kegg_aliases",1],
			["metacyc_aliases",1],
			["names",1],
			["ec_numbers",1],
			["roles",2,";"]
		];
		for (my $k=0; $k < @{$headings}; $k++) {
			if (defined($rxnobj->{$headings->[$k]->[0]})) {
				if ($headings->[$k]->[1] == 0) {
					$line .= $rxnobj->{$headings->[$k]->[0]};
				} elsif ($headings->[$k]->[1] == 1) {
					$line .= join("|",@{$rxnobj->{$headings->[$k]->[0]}});
				} else {
					my $subarray = [];
					for (my $m=0; $m < @{$rxnobj->{$headings->[$k]->[0]}}; $m++) {
						my $miniarray = [split(/;/,$rxnobj->{$headings->[$k]->[0]}->[$m])];
						push(@{$subarray},$miniarray->[1]);
					}
					$line .= join("|",@{$subarray});
				}
			}
			$line .= "\t";
		}
	} else {
		$line .= "\t\t\t\t\t\t\t\t\t";
	}
	#Now printing strain specific data
	my $priority = 0;
	my $rxnobj = $rxnmodeldata->{$id};
	for (my $i=0; $i < @{$genomelist}; $i++) {
		if ($i > 0) {
			$line .= "\t";
		}
		if (defined($rxnobj->{$genomelist->[$i]})) {
			my $rxngdata = $rxnobj->{$genomelist->[$i]};
			if (length($rxngdata->{gapfillString}) > 0) {
				$line .= "G/";
			}
			if ($rxngdata->{ssflux} > 0) {
				$line .= "PosF/";
			} elsif ($rxngdata->{ssflux} < 0) {
				$line .= "NegF/";
			} else{
				$line .= "NoF/";
			}
			$line .= $rxngdata->{mmclass}."/".$rxngdata->{comclass}."/".$rxngdata->{direction}."/";
			if (defined($rxngdata->{biomassdep})) {
				for (my $j=0; $j < @{$rxngdata->{biomassdep}};  $j++) {
					if ($j > 0) {
						$line .= ";";
					}
					my $cpdid = $rxngdata->{biomassdep}->[$j]->[1];
					$cpdid =~ s/_c0//;
					$line .= $globalcpdhash->{$cpdid}->{name};
				}
				$line .= "/";
			}
			my $genelists = "";
			if (defined($rxngdata->{genelist})) {
				$genelists = join(",",@{$rxngdata->{genelist}});
			}
			$genelists =~ s/fig\|\d+\.\d+\.peg\.//g;
			$line .= $genelists."/";
			if ($rxngdata->{ssflux} != 0) {
				$priority += 10;
			} elsif ($rxngdata->{mmclass} eq "p" && $rxngdata->{mmclass} eq "n") {
				$priority += 8;
			} elsif ($rxngdata->{mmclass} ne "b") {
				$priority += 6;
			} elsif ($rxngdata->{comclass} ne "b") {
				$priority += 4;
			} elsif (length($rxngdata->{gapfillString}) > 0) {
				$priority += 2;
			}
		} else {
			$line .= "-";
		}
	}
	return ($line,$priority);
}
