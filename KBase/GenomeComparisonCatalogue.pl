#!/usr/bin/perl -w

use strict;
use JSON::XS;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
$|=1;

my $directory = $ARGV[0];
(my $protcomp,my $meta) = get_workspace_object("chenry:1452494413455/MGH7857_KPPR1_comparison");
(my $kppr1,$meta) = get_workspace_object("chenry:1452494413455/Klebsiella_pneumoniae_kppr1_kbase");
(my $mgh,$meta) = get_workspace_object("chenry:1452494413455/Klebsiella_pneumoniae_MGH78578_kbase");
(my $pangenome,$meta) = get_workspace_object("chenry:1452494413455/MGH78578_KPPR1_protein_families");

my $ftrkppr1 = {};
my $ftrmgh = {};

my $functions = {};

my $mgh_orphan_roles = {};
my $kppr1_orphan_roles = {};

for (my $i=0; $i < @{$kppr1->{features}}; $i++) {
	$ftrkppr1->{$kppr1->{features}->[$i]->{id}} = $kppr1->{features}->[$i];
	my $roles = [split(/\s*;\s+|\s+[\@\/]\s+/,$kppr1->{features}->[$i]->{function})];
	for (my $j=0; $j < @{$roles}; $j++) {
		$functions->{$roles->[$j]}->{kpprftrs}->{$kppr1->{features}->[$i]->{id}} = 1;
	}
}
for (my $i=0; $i < @{$mgh->{features}}; $i++) {
	$ftrmgh->{$mgh->{features}->[$i]->{id}} = $mgh->{features}->[$i];
	my $roles = [split(/\s*;\s+|\s+[\@\/]\s+/,$mgh->{features}->[$i]->{function})];
	for (my $j=0; $j < @{$roles}; $j++) {
		$functions->{$roles->[$j]}->{mghftrs}->{$mgh->{features}->[$i]->{id}} = 1;
	}
}

my $data = $protcomp->{data1};
my $list = $protcomp->{proteome1names};
my $olist = $protcomp->{proteome2names};
my $translate;
for(my $i=0; $i < @{$data}; $i++) {
	for (my $j=0; $j < @{$data->[$i]}; $j++) {
		if ($data->[$i]->[$j]->[2] == 100) {
			if (defined($ftrkppr1->{$list->[$i]})) {
				$ftrkppr1->{$list->[$i]}->{orthos}->{$olist->[$data->[$i]->[$j]->[0]]} = "blast";
				$ftrmgh->{$olist->[$data->[$i]->[$j]->[0]]}->{orthos}->{$list->[$i]} = "blast";
			} elsif (defined($ftrmgh->{$list->[$i]})) {
				$ftrmgh->{$list->[$i]}->{orthos}->{$olist->[$data->[$i]->[$j]->[0]]} = "blast";
				$ftrkppr1->{$olist->[$data->[$i]->[$j]->[0]]}->{orthos}->{$list->[$i]} = "blast";
			}
		}
	}
}

for(my $i=0; $i < @{$pangenome->{orthologs}}; $i++) {
	my $ortho = $pangenome->{orthologs}->[$i];
	for(my $j=0; $j < @{$ortho->{orthologs}}; $j++) {
		for(my $k=0; $k < @{$ortho->{orthologs}}; $k++) {
			if ($j != $k) {
				if (defined($ftrkppr1->{$ortho->{orthologs}->[$j]->[0]})) {
					if (defined($ftrmgh->{$ortho->{orthologs}->[$k]->[0]})) {
						if (defined($ftrkppr1->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]})) {
							$ftrkppr1->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]} .= ";kmer";
						} else {
							$ftrkppr1->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]} = "kmer";
						}
					} elsif (defined($ftrkppr1->{$ortho->{orthologs}->[$k]->[0]})) {
						$ftrkppr1->{$ortho->{orthologs}->[$j]->[0]}->{paras}->{$ortho->{orthologs}->[$k]->[0]} = "kmer";
					}
				} elsif (defined($ftrmgh->{$ortho->{orthologs}->[$j]->[0]})) {
					if (defined($ftrkppr1->{$ortho->{orthologs}->[$k]->[0]})) {
						if (defined($ftrmgh->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]})) {
							$ftrmgh->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]} .= ";kmer";
						} else {
							$ftrmgh->{$ortho->{orthologs}->[$j]->[0]}->{orthos}->{$ortho->{orthologs}->[$k]->[0]} = "kmer";
						}
					} elsif (defined($ftrmgh->{$ortho->{orthologs}->[$k]->[0]})) {
						$ftrmgh->{$ortho->{orthologs}->[$j]->[0]}->{paras}->{$ortho->{orthologs}->[$k]->[0]} = "kmer";
					}
				}
			}
		}
	}
}

foreach my $ftr (keys(%{$ftrkppr1})) {
	if (!defined($ftrkppr1->{$ftr}->{orthos})) {
		my $roles = [split(/\s*;\s+|\s+[\@\/]\s+/,$ftrkppr1->{$ftr}->{function})];
		for (my $j=0; $j < @{$roles}; $j++) {
			$functions->{$roles->[$j]}->{kppr1_noorthos}->{$ftr} = 1;
		}
	}
}
foreach my $ftr (keys(%{$ftrmgh})) {
	if (!defined($ftrmgh->{$ftr}->{orthos})) {
		my $roles = [split(/\s*;\s+|\s+[\@\/]\s+/,$ftrmgh->{$ftr}->{function})];
		for (my $j=0; $j < @{$roles}; $j++) {
			$functions->{$roles->[$j]}->{mgh_noorthos}->{$ftr} = 1;
		}
	}
}


my $filename = "/Users/chenry/code/ModelSEEDDatabase/SOLRDump/Subsystems.json";
open(my $fhh, "<", $filename);
my $jsonfile = "";
while (my $line = <$fhh>) {
	chomp($line);
	$jsonfile .= $line;
}
my $subsysdata = decode_json $jsonfile;
close($fhh);
my $ssdatahash;
foreach my $ss (@{$subsysdata}) {
	$ssdatahash->{$ss->{id}} = $ss;
}
$filename = "/Users/chenry/code/ModelSEEDDatabase/SOLRDump/Roles.json";
open(my $fh, "<", $filename);
$jsonfile = "";
while (my $line = <$fh>) {
	chomp($line);
	$jsonfile .= $line;
}
close($fh);
my $roledata = decode_json $jsonfile;
my $rolehash;
foreach my $role (@{$roledata}) {
	if ($role->{source} ne "KEGG") {
		$rolehash->{$role->{searchname}} = $role;
	}
}

my $outfile = $directory."Functions.tsv";
open(my $fo, ">", $outfile);
print $fo "Function\tSubsystem\tClass\tReactions\tMGH genes\tKPPR1 genes\tMGH orthos\tKPPR1 orthos\n";
my $sshash = {};
my $classhash = {};
foreach my $func (keys(%{$functions})) {
	my $searchname = $func;
	$searchname = lc($searchname);
	$searchname =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
	$searchname =~ s/\s//g;
	$searchname =~ s/\#.*$//g;
	$searchname =~ s/\(ec\)//g;
	print $fo $func."\t";
	if (defined($rolehash->{$searchname})) {
		my $funcss;
		my $funcclasses;
		my $funcrxns;
		foreach my $rxn (@{$rolehash->{$searchname}->{reactions}}) {
			my $array = [split(/;/,$rxn)];
			$funcrxns->{$array->[0]} = 1;
		}
		foreach my $ss (@{$rolehash->{$searchname}->{subsystems}}) {
			my $array = [split(/;/,$ss)];
			if ($array->[2] ne "Experimental Subsystems" && $array->[2] ne "Clustering-based subsystems") {
				$funcss->{$array->[1]} = 1;
				$funcclasses->{$array->[2]} = 1;
				if (!defined($sshash->{$array->[0]})) {
					$sshash->{$array->[0]} = {
						id => $array->[0],
						name => $array->[1],
						class => $array->[2],
						subclass => $array->[3],
						function_count => 0,
						shared => 0,
						kppr_only => 0,
						mgh_only => 0,
						kppr_no_orthos => 0,
						mgh_no_orthos => 0,
						metaroles => 0
					};
					if (defined($ssdatahash->{$array->[0]})) {
						$sshash->{$array->[0]}->{function_count} = @{$ssdatahash->{$array->[0]}->{roles}};
					}
				}
				if (!defined($classhash->{$array->[2]})) {
					$classhash->{$array->[2]} = {
						name => $array->[2],
						subsyshash => {},
						function_count => 0,
						shared => 0,
						kppr_only => 0,
						mgh_only => 0,
						kppr_no_orthos => 0,
						mgh_no_orthos => 0,
						metaroles => 0
					};
				}
				if (!defined($classhash->{$array->[2]}->{subsyshash}->{$array->[0]})) {
					$classhash->{$array->[2]}->{subsyshash}->{$array->[0]} = 1;
					$classhash->{$array->[2]}->{function_count} += $sshash->{$array->[0]}->{function_count};
				}
				if (defined($functions->{$func}->{mghftrs}) && defined($functions->{$func}->{kpprftrs})) {
					$classhash->{$array->[2]}->{shared}++;
					$sshash->{$array->[0]}->{shared}++;
				} elsif (defined($functions->{$func}->{mghftrs})) {
					$classhash->{$array->[2]}->{mgh_only}++;
					$sshash->{$array->[0]}->{mgh_only}++;
				} elsif (defined($functions->{$func}->{kpprftrs})) {
					$classhash->{$array->[2]}->{kppr_only}++;
					$sshash->{$array->[0]}->{kppr_only}++;
				}
				if (defined($functions->{$func}->{kppr1_noorthos})) {
					$classhash->{$array->[2]}->{kppr_no_orthos}++;
					$sshash->{$array->[0]}->{kppr_no_orthos}++;
				}
				if (defined($functions->{$func}->{mgh_noorthos})) {
					$classhash->{$array->[2]}->{mgh_no_orthos}++;
					$sshash->{$array->[0]}->{mgh_no_orthos}++;
				}
				if (keys(%{$funcrxns}) > 0) {
					$classhash->{$array->[2]}->{metaroles}++;
					$sshash->{$array->[0]}->{metaroles}++;
				}
			}
		}
		print $fo join("|",keys(%{$funcss}))."\t".join("|",keys(%{$funcclasses}))."\t".join("|",keys(%{$funcrxns}))."\t";
	} else {
		print $fo "\t\t\t";
	}
	if (defined($functions->{$func}->{mghftrs})) {
		print $fo join(";",keys(%{$functions->{$func}->{mghftrs}}));
	}
	print $fo "\t";
	if (defined($functions->{$func}->{kpprftrs})) {
		print $fo join(";",keys(%{$functions->{$func}->{kpprftrs}}));
	}
	print $fo "\t";
	if (defined($functions->{$func}->{mgh_noorthos})) {
		print $fo join(";",keys(%{$functions->{$func}->{mgh_noorthos}}));
	}
	print $fo "\t";
	if (defined($functions->{$func}->{kppr1_noorthos})) {
		print $fo join(";",keys(%{$functions->{$func}->{kppr1_noorthos}}));
	}
	print $fo "\n";
}
close($fo);

$outfile = $directory."Subsystems.tsv";
open(my $foo, ">", $outfile);
print $foo "ID\tSubsystem\tClass\tSubclass\tFunction count\tMetabolic roles\tShared functions\tKPPR1 genes only\tMGH genes only\tKPPR1 no orthos\tMGH no orthos\n";
foreach my $ss (keys(%{$sshash})) {
	if ($sshash->{$ss}->{class} ne "Experimental Subsystems" && $sshash->{$ss}->{class} ne "Clustering-based subsystems") {
		print $foo $sshash->{$ss}->{id}."\t".$sshash->{$ss}->{name}."\t".$sshash->{$ss}->{class}."\t".$sshash->{$ss}->{subclass}."\t".$sshash->{$ss}->{function_count}."\t".$sshash->{$ss}->{metaroles}."\t".$sshash->{$ss}->{shared}."\t".$sshash->{$ss}->{kppr_only}."\t".$sshash->{$ss}->{mgh_only}."\t".$sshash->{$ss}->{kppr_no_orthos}."\t".$sshash->{$ss}->{mgh_no_orthos}."\n";
	}
}
close($foo);

$outfile = $directory."SubsystemClass.tsv";
open(my $fooo, ">", $outfile);
print $fooo "Class\tFunction count\tMetabolic roles\tShared functions\tKPPR1 genes only\tMGH genes only\tMGH no orthos\tKPPR1 no orthos\n";
foreach my $class (keys(%{$classhash})) {
	print $fooo $classhash->{$class}->{name}."\t".$classhash->{$class}->{function_count}."\t".$classhash->{$class}->{metaroles}."\t".$classhash->{$class}->{shared}."\t".$classhash->{$class}->{kppr_only}."\t".$classhash->{$class}->{mgh_only}."\t".$classhash->{$class}->{kppr_no_orthos}."\t".$classhash->{$class}->{mgh_no_orthos}."\n";
}
close($fooo);