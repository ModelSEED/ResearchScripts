#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

my $cpxs = {};
my $rolehash = {};
my $rxnprophash = {};
my $transhash = {};
open(my $ffh, "<", $directory."ReactionData.tsv");
my $line = <$ffh>
while ($line = <$ffh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if ($array->[5] == 1) {
		$transhash->{$array->[0]} = 1;
	}
}
close($ffh);

open(my $fh, "<", $directory."Roles.tsv");
$line = <$fh>
while ($line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$array->[1] = lc($array->[1]);
	$array->[1] =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
	$array->[1] =~ s/\s//g;
	$array->[1] =~ s/\#.*$//g;
	$array->[1] =~ s/\(ec\)//g;
	$rolehash->{$array->[0]} = $array->[1];
}
close($fh);

open(my $fhh, "<", $directory."Complexes.tsv");
$line = <$fhh>
while ($line = <$fhh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	my $roles = [split(/\|/,$array->[5])];
	for (my $i=0; $i < @{$roles}; $i++) {
		my $items = [split(/;/,$roles->[$i])];
		if (defined($rolehash->{$items->[0]})) {
			$cpxs->{$array->[0]}->{$rolehash->{$items->[0]}} = 1;
		}
	}
}
close($fhh);

my $rxnhash;
my $count = 0;
open(my $fhhh, "<", $directory."TemplateReactions.tsv");
$line = <$fhhh>
while ($line = <$fhhh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	my $cpxsarray = [split(/\|/,$array->[8])];
	if (!defined($rxnprophash->{$array->[0]})) {
		$rxnprophash->{$array->[0]} = {
			cpx => 0,
			transport => 0,
			models => {modelseed => 1}
		};
		if (defined($transhash->{$array->[0]})) {
			$rxnprophash->{$array->[0]}->{transport} = 1;
		}
	}
	for (my $i=0; $i < @{$cpxsarray}; $i++) {
		if (defined($cpxs->{$cpxsarray->[$i]})) {
			my $rolecount = keys(%{$cpxs->{$cpxsarray->[$i]}});
			if ($rolecount > 1) {
				$rxnprophash->{$array->[0]}->{cpx} = 1;
			}
			foreach my $rolename (keys(%{$cpxs->{$cpxsarray->[$i]}})) {
				if (!defined($rxnhash->{$array->[0]}->{$rolename})) {
					$rxnhash->{$array->[0]}->{$rolename} = 0;
				}
				$rxnhash->{$array->[0]}->{$rolename} += 10000;
				$count++;
			}
		}
	}
}
close($fhhh);
print "COUNT:".$count."\n";
my $genomes;
($genomes->{iJO1366},my $info) = get_workspace_object("10559/12/1");
($genomes->{iRsp1140},$info) = get_workspace_object("10559/11/1");
($genomes->{iMM904},$info) = get_workspace_object("10559/16/1");
($genomes->{iJN678},$info) = get_workspace_object("10559/10/1");
($genomes->{iYL1228},$info) = get_workspace_object("10559/13/1");
($genomes->{iBsu1103},$info) = get_workspace_object("10559/15/1");
($genomes->{iAbaylyiV4},$info) = get_workspace_object("10559/17/1");
($genomes->{iAH991},$info) = get_workspace_object("10559/14/1");
my $genehash;
foreach my $model (keys(%{$genomes})) {
	my $genome = $genomes->{$model};
	my $ftrs = $genome->{features};
	for (my $i = 0; $i < @{$ftrs}; $i++) {
		my $function = $ftrs->[$i]->{function};
		if (defined($function)) {
			$function =~ s/\s*#.+//;
			my $fnarray = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$function))];
			for (my $j=0; $j < @{$fnarray}; $j++) {
				$fnarray->[$j] = lc($fnarray->[$j]);
				$fnarray->[$j] =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
				$fnarray->[$j] =~ s/\s//g;
				$fnarray->[$j] =~ s/\#.*$//g;
				$fnarray->[$j] =~ s/\(ec\)//g;
				$genehash->{$ftrs->[$i]->{id}}->{$fnarray->[$j]} = 1;
			}
		} else {
			print "No function!\n";
		}
	}
}

$genomes = {};
($genomes->{iJO1366},$info) = get_workspace_object("8380/3/1");
($genomes->{iRsp1140},$info) = get_workspace_object("8380/1/2");
($genomes->{iMM904},$info) = get_workspace_object("8380/6/4");
($genomes->{iJN678},$info) = get_workspace_object("8380/2/5");
($genomes->{iYL1228},$info) = get_workspace_object("8380/8/7");
($genomes->{iBsu1103},$info) = get_workspace_object("8380/4/4");
($genomes->{iAbaylyiV4},$info) = get_workspace_object("8380/5/4");
($genomes->{iAH991},$info) = get_workspace_object("8380/7/11");
foreach my $model (keys(%{$genomes})) {
	my $genome = $genomes->{$model};
	my $rxns = $genome->{modelreactions};
	for (my $i=0;$i < @{$rxns}; $i++) {
		my $id = $rxns->[$i]->{id};
		if ($rxns->[$i]->{reaction_ref} =~ m/\/([^\/]+)$/) {
			if ($1 ne "rxn00000") {
				$id = $1;
			}
		}
		if ($id =~ m/(.+)_[a-z]\d+$/) {
			$id = $1;
		}
		if (!defined($rxnprophash->{$id})) {
			$rxnprophash->{$id} = {
				cpx => 0,
				transport => 0,
				models => {}
			};
		}
		$rxnprophash->{$id}->{models}->{$model} = 1;
		if (defined($transhash->{$id})) {
			$rxnprophash->{$id}->{transport} = 1;
		} else {
			my $rgts = $rxns->[$i]->{modelReactionReagents};
			my $comphash = {};
			for (my $j=0;$j < @{$rgts}; $j++) {
				if ($rgts->[$j]->{modelcompound_ref} =~ m/_([a-z]\d+)$/) {
					$comphash->{$1} = 1;
				}
			}
			if (keys(%{$comphash}) > 1) {
				$rxnprophash->{$id}->{transport} = 1;
			}
		}
		my $prots = $rxns->[$i]->{modelReactionProteins};
		for (my $j=0;$j < @{$prots}; $j++) {
			my $subunits = $prots->[$j]->{modelReactionProteinSubunits};
			if (@{$subunits} > 1) {
				$rxnprophash->{$id}->{cpx} = 1;
			}
			for (my $k=0;$k < @{$subunits}; $k++) {
				my $features = $subunits->[$k]->{feature_refs};
				for (my $m=0;$m < @{$features}; $m++) {
					if ($features->[$m] =~ m/\/([^\/]+)$/) {
						my $geneid = $1;
						if (defined($genehash->{$geneid})) {
							foreach my $func (keys(%{$genehash->{$geneid}})) {
								if (!defined($rxnhash->{$id}->{$func})) {
									$rxnhash->{$id}->{$func} = 0;
								}
								$rxnhash->{$id}->{$func}++;
							}
						}
					}
				}
			}
		}
	}
}
open (my $out, ">", $directory."/ReactionFunctions.txt");
print $out "Reaction\tFunction\tCount\tComplex\tTransport\tModels\n";
foreach my $key (keys(%{$rxnhash})) {
	foreach my $func (keys(%{$rxnhash->{$key}})) {
		print $out $key."\t".$func."\t".$rxnhash->{$key}->{$func}."\t".$rxnprophash->{$key}->{cpx}."\t".$rxnprophash->{$key}->{transport}."\t".join("|",keys(%{$rxnprophash->{$key}->{models}}))."\n";
	}
}
close($out);