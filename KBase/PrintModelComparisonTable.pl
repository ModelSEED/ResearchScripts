#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use POSIX qw/floor/;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $ws = Bio::KBase::kbaseenv::ws_client();

my $directory = $ARGV[0];
my $workspace = $ARGV[1];

my $genomes = [qw(
100226.1
192222.1
221988.1
243265.1
266940.5
290402.34
326442.4
62977.3
106370.11
196600.1
223926.1
243273.1
267608.1
290633.1
335992.3
71421.1
107806.1
196627.4
224308.1
243274.1
267671.1
292415.3
342108.5
76114.4
115711.7
198094.1
224324.1
243276.1
292459.1
345073.6
76869.3
122586.1
198215.1
224326.1
243277.1
269798.12
295405.3
360095.3
83332.1
158878.1
203267.1
224911.1
246194.3
269799.3
296591.1
360110.3
83333.1
158879.1
203907.1
224914.1
246200.3
272560.3
298386.1
360111.3
85962.1
160490.1
205922.3
227377.1
247156.1
272561.1
299768.3
36870.1
85963.1
160492.1
206672.1
228410.1
257313.1
272562.1
300852.3
36873.1
93061.3
163164.1
207559.3
234826.3
258594.1
272623.1
302409.3
370552.3
99287.1
169963.1
208964.1
235909.3
261594.1
272624.3
309807.19
376686.6
170187.1
211586.1
262768.1
272626.1
316407.3
386656.4
171101.1
212717.1
242619.1
264203.3
272635.1
393117.3
176299.3
214092.1
243164.3
264462.1
272947.1
319701.3
393130.3
177416.3
216596.1
243230.1
265072.7
281090.3
323098.3
399726.4
190304.1
216598.1
243231.1
266117.6
283942.3
323261.3
401614.5
190650.1
220668.1
243233.4
266834.1
290397.13
326298.3
445932.3
)];

my $letters = [qw(
A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
)];

my $rxnhash = Bio::KBase::utilities::reaction_hash();
my $genehash;
my $geneanno;
my $unique_combinations_hash;
#open (my $fout, ">", $directory."/FullTable.txt");
#print $fout "Genome\tGene\tRxn\tEquation\tEC numbers\tCurrent roles\tOrig\tRAST\tRAST2\tRAST annotation\tRAST2 annotation\n";
foreach my $genome (@{$genomes}) {
	eval {
	my $g1 = $handler->util_get_object($workspace."/".$genome.".RAST",{raw => 1});
	my $g2 = $handler->util_get_object($workspace."/".$genome.".RAST2",{raw => 1});
	my $m1 = $handler->util_get_object($workspace."/Seed".$genome,{raw => 1});
	my $m2 = $handler->util_get_object($workspace."/".$genome.".RAST.fbamodel",{raw => 1});
	my $m3 = $handler->util_get_object($workspace."/".$genome."RAST2.fbamodel",{raw => 1});
	my $genehash = {};
	my $geneanno = {};
	for (my $i=0; $i < @{$g1->{features}}; $i++) {
		if ($g1->{features}->[$i]->{id} =~ m/(peg\.\d+)/) {
			my $geneid = $1;
			$geneanno->{$geneid}->{RAST} = $g1->{features}->[$i]->{function};
		}
	}
	for (my $i=0; $i < @{$g2->{features}}; $i++) {
		if ($g2->{features}->[$i]->{id} =~ m/(peg\.\d+)/) {
			my $geneid = $1;
			$geneanno->{$geneid}->{RAST2} = $g2->{features}->[$i]->{function};
		}
	}
	for (my $i=0; $i < @{$m1->{modelreactions}}; $i++) {
		if ($m1->{modelreactions}->[$i]->{id} =~ m/(rxn\d+)/) {
			my $rxnid = $1;
			my $prots = $m1->{modelreactions}->[$i]->{modelReactionProteins};
			for (my $j=0; $j < @{$prots}; $j++) {
				my $subunits = $prots->[$j]->{modelReactionProteinSubunits};
				for (my $k=0; $k < @{$subunits}; $k++) {
					for (my $m=0; $m < @{$subunits->[$k]->{feature_refs}}; $m++) {
						if ($subunits->[$k]->{feature_refs}->[$m] =~ m/(peg\.\d+)/) {
							my $geneid = $1;
							$genehash->{$geneid}->{$rxnid}->{orig} = 1;
						}
					}
				}
			}
		}
	}
	for (my $i=0; $i < @{$m2->{modelreactions}}; $i++) {
		if ($m2->{modelreactions}->[$i]->{id} =~ m/(rxn\d+)/) {
			my $rxnid = $1;
			my $prots = $m2->{modelreactions}->[$i]->{modelReactionProteins};
			for (my $j=0; $j < @{$prots}; $j++) {
				my $subunits = $prots->[$j]->{modelReactionProteinSubunits};
				for (my $k=0; $k < @{$subunits}; $k++) {
					for (my $m=0; $m < @{$subunits->[$k]->{feature_refs}}; $m++) {
						if ($subunits->[$k]->{feature_refs}->[$m] =~ m/(peg\.\d+)/) {
							my $geneid = $1;
							$genehash->{$geneid}->{$rxnid}->{RAST} = 1;
						}
					}
				}
			}
		}
	}
	for (my $i=0; $i < @{$m3->{modelreactions}}; $i++) {
		if ($m3->{modelreactions}->[$i]->{id} =~ m/(rxn\d+)/) {
			my $rxnid = $1;
			my $prots = $m3->{modelreactions}->[$i]->{modelReactionProteins};
			for (my $j=0; $j < @{$prots}; $j++) {
				my $subunits = $prots->[$j]->{modelReactionProteinSubunits};
				for (my $k=0; $k < @{$subunits}; $k++) {
					for (my $m=0; $m < @{$subunits->[$k]->{feature_refs}}; $m++) {
						if ($subunits->[$k]->{feature_refs}->[$m] =~ m/(peg\.\d+)/) {
							my $geneid = $1;
							$genehash->{$geneid}->{$rxnid}->{RAST2} = 1;
						}
					}
				}
			}
		}
	}
	foreach my $gene (keys(%{$genehash})) {
		foreach my $rxnid (keys(%{$genehash->{$gene}})) {
			my $present = [0,0,0];
			if (defined($genehash->{$gene}->{$rxnid}->{orig})) {
				$present->[0] = 1;
			}
			if (defined($genehash->{$gene}->{$rxnid}->{RAST})) {
				$present->[1] = 1;
			}
			if (defined($genehash->{$gene}->{$rxnid}->{RAST2})) {
				$present->[2] = 1;
			}
			$unique_combinations_hash->{$rxnid}->{$present->[0]}->{$present->[1]}->{$present->[2]}->{$geneanno->{$gene}->{RAST}}->{$geneanno->{$gene}->{RAST2}}->{$genome}->{$gene} = 1;
			my $eqn = "";
			my $ec = "";
			my $roles = "";
			if (defined($rxnhash->{$rxnid})) {
				$eqn = $rxnhash->{$rxnid}->{definition};
				if (defined($rxnhash->{$rxnid}->{ec_numbers})) {
					$ec = join("|",@{$rxnhash->{$rxnid}->{ec_numbers}});
				}
				if (defined($rxnhash->{$rxnid}->{roles})) {
					for (my $i=0; $i < @{$rxnhash->{$rxnid}->{roles}}; $i++) {
						my $temparray = [split/;/,$rxnhash->{$rxnid}->{roles}->[$i]];
						if (length($roles) > 0) {
							$roles .= "|";
						}
						$roles .= $temparray->[1];
					}
				}
			}
			if ($present->[0] == 0 || $present->[1] == 0 || $present->[2] == 0) {
				#print $fout $genome."\t".$gene."\t".$rxnid."\t".$eqn."\t".$ec."\t".$roles."\t";
				#print $fout $present->[0]."\t".$present->[1]."\t".$present->[2]."\t".$geneanno->{$gene}->{RAST}."\t".$geneanno->{$gene}->{RAST2}."\n";
			}
		}
	}
	};
	if ($@) {
		print "Error on ".$genome."\n".$@."\n";
	}
}
#close($fout);
my $alllabels = {};
open (my $fout2, ">", $directory."/CombinedTable.txt");
print $fout2 "Num genomes\tNum genes\tRxn\tEquation\tEC numbers\t111\t011\t001\t010\t100\t110\t101\tCurrent roles\tAnnotations\n";
foreach my $rxnid (keys(%{$unique_combinations_hash})) {
	my $instances;
	my $rolehash = {};
	my $totalgenome = 0;
	my $totalgenes = 0;
	foreach my $orig (keys(%{$unique_combinations_hash->{$rxnid}})) {
		foreach my $rast (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}})) {
			foreach my $rast2 (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}->{$rast}})) {
				$genehash = {};
				my $genomehash = {};
				foreach my $rastrole (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}->{$rast}->{$rast2}})) {
					if (!defined($rolehash->{$rastrole}->{newcount})) {
						$rolehash->{$rastrole}->{count} = 0;
						$rolehash->{$rastrole}->{newcount} = {
							rast1 => 0,
							all => 0,
							rast2 => 0
						};
					}
					foreach my $rast2role (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}->{$rast}->{$rast2}->{$rastrole}})) {
						if (!defined($rolehash->{$rast2role}->{newcount})) {
							$rolehash->{$rast2role}->{count} = 0;
							$rolehash->{$rast2role}->{newcount} = {
								rast1 => 0,
								all => 0,
								rast2 => 0
							};
						}
						foreach my $genome (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}->{$rast}->{$rast2}->{$rastrole}->{$rast2role}})) {
							$genomehash->{$genome} = 1; 
							foreach my $gene (keys(%{$unique_combinations_hash->{$rxnid}->{$orig}->{$rast}->{$rast2}->{$rastrole}->{$rast2role}->{$genome}})) {
								$genehash->{$gene} = 1;
								if (!defined($rolehash->{$rastrole}->{$orig.$rast.$rast2})) {
									$rolehash->{$rastrole}->{$orig.$rast.$rast2} = {
										rast1 => 0,
										all => 0,
										rast2 => 0
									};
								}
								$rolehash->{$rastrole}->{$orig.$rast.$rast2}->{all}++;
								$rolehash->{$rastrole}->{$orig.$rast.$rast2}->{rast1}++;
								if ($orig.$rast.$rast2 eq "001" || $orig.$rast.$rast2 eq "100" || $orig.$rast.$rast2 eq "101") {
									$rolehash->{$rastrole}->{newcount}->{all}++;
									$rolehash->{$rastrole}->{newcount}->{rast1}++;
								}
								if (!defined($rolehash->{$rast2role}->{$orig.$rast.$rast2})) {
									$rolehash->{$rast2role}->{$orig.$rast.$rast2} = {
										rast1 => 0,
										all => 0,
										rast2 => 0
									};
								}
								$rolehash->{$rastrole}->{count}++;
								$rolehash->{$rast2role}->{count}++;
								$rolehash->{$rast2role}->{$orig.$rast.$rast2}->{all}++;
								$rolehash->{$rast2role}->{$orig.$rast.$rast2}->{rast2}++;
								if ($orig.$rast.$rast2 eq "010" || $orig.$rast.$rast2 eq "100" || $orig.$rast.$rast2 eq "110") {
									$rolehash->{$rast2role}->{newcount}->{all}++;
									$rolehash->{$rast2role}->{newcount}->{rast2}++;
								}
							}
						}
					}
				}	
				my $genecount = keys(%{$genehash});
				my $genomecount = keys(%{$genomehash});
				$instances->{$orig.$rast.$rast2} = [$genomecount,$genecount];
				$totalgenome += $genomecount;
				$totalgenes += $genecount;
			}
		}
	}
	my $types = ["111","011","001","010","100","110","101"];
	my $eqn = "";
	my $ec = "";
	my $roles = "";
	if (defined($rxnhash->{$rxnid})) {
		$eqn = $rxnhash->{$rxnid}->{definition};
		if (defined($rxnhash->{$rxnid}->{ec_numbers})) {
			$ec = join("|",@{$rxnhash->{$rxnid}->{ec_numbers}});
		}
		if (defined($rxnhash->{$rxnid}->{roles})) {
			for (my $i=0; $i < @{$rxnhash->{$rxnid}->{roles}}; $i++) {
				my $temparray = [split/;/,$rxnhash->{$rxnid}->{roles}->[$i]];
				if (length($roles) > 0) {
					$roles .= "|";
				}
				my $count = 0;
				if (defined($rolehash->{$temparray->[1]})) {
					$count = $rolehash->{$temparray->[1]}->{count};
				}
				my $id = $i+1;
				$alllabels->{$rxnid}->{currroles}->{$temparray->[1]} = $id;
				$roles .= $id.":".$temparray->[1].":".$count;
			}
		}
	}
	print $fout2 $totalgenome."\t".$totalgenes."\t".$rxnid."\t".$eqn."\t".$ec;
	foreach my $type (@{$types}) {
		if (defined($instances->{$type})) {
			print $fout2 "\t".$instances->{$type}->[1];
		} else {
			print $fout2 "\t0";
		}
	}
	print $fout2 "\t".$roles;
	my $rolelist = [sort { $rolehash->{$b}->{newcount}->{all} <=> $rolehash->{$a}->{newcount}->{all} } keys(%{$rolehash})];
	my $count = 0;
	foreach my $role (@{$rolelist}) {
		if ($rolehash->{$role}->{newcount}->{all} > 0) {
			my $lettercount = floor($count/26)+1;
			my $mod = $count % 26;
			my $label = "";
			for (my $i=0; $i < $lettercount; $i++) {
				$label .= $letters->[$mod];
			}
			$count++;
			$alllabels->{$rxnid}->{newroles}->{$role} = $label;
			print $fout2 "\t".$label."|".$rolehash->{$role}->{newcount}->{rast1}."|".$rolehash->{$role}->{newcount}->{rast2}.";".$role.";";
			if (defined($rolehash->{$role}->{"001"})) {
				print $fout2 $rolehash->{$role}->{"001"}->{rast1}.":";
			} else {
				print $fout2  "0:";
			}
			if (defined($rolehash->{$role}->{"100"})) {
				print $fout2 $rolehash->{$role}->{"100"}->{rast1}.":";
			} else {
				print $fout2  "0:";
			}
			if (defined($rolehash->{$role}->{"101"})) {
				print $fout2 $rolehash->{$role}->{"101"}->{rast1};
			} else {
				print $fout2  "0";
			}
			print $fout2 "|";
			if (defined($rolehash->{$role}->{"010"})) {
				print $fout2 $rolehash->{$role}->{"010"}->{rast2}.":";
			} else {
				print $fout2  "0:";
			}
			if (defined($rolehash->{$role}->{"100"})) {
				print $fout2 $rolehash->{$role}->{"100"}->{rast2}.":";
			} else {
				print $fout2  "0:";
			}
			if (defined($rolehash->{$role}->{"110"})) {
				print $fout2 $rolehash->{$role}->{"110"}->{rast2};
			} else {
				print $fout2  "0";
			}
		}
	}
	print $fout2 "\n";	
}
close($fout2);
Bio::KBase::ObjectAPI::utilities::PRINTFILE($directory."/RoleTables.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($alllabels,1)]);