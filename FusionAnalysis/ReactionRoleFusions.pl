#!/usr/bin/perl

use strict;
use warnings;

my $directory = $ARGV[0];

my $rxnacthash = {};
open(my $fh, "<", $directory."/ReactionActivityAnalysis.txt");
my $line = <$fh>;
while ($line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$rxnacthash->{$items->[0]."\t".$items->[1]} = $items;
}
close($fh);

my $functionhash;
open(my $fffh, "<", $directory."/FunctionFusionAnalysisTwo.txt");
$line = <$fffh>;
while ($line = <$fffh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$functionhash->{$items->[0]} = {
		genes => $items->[1],
		genomes => $items->[2],
		fusions => $items->[3],
		fusiongenomes => $items->[4],
		reactions => {}
	};
}
close($fffh);

my $reactionhash;
open(my $fhh, "<", $directory."/ReactionFunctions.txt");
$line = <$fhh>;
while ($line = <$fhh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	if (!defined($reactionhash->{$items->[0]})) {
		$reactionhash->{$items->[0]} = {
			cpx => $items->[3],
			transport => $items->[4],
			models => $items->[5],
			functions => {}
		};
	}
	if (!defined($reactionhash->{$items->[0]}->{functions}->{$items->[1]})) {
		$reactionhash->{$items->[0]}->{functions}->{$items->[1]} = 0;
		if (!defined($functionhash->{$items->[1]})) {
			$functionhash->{$items->[1]} = {
				genes => 1,
				genomes => 1,
				fusions => 0,
				fusiongenomes => 0,
				reactions => {}
			};
		}
		$functionhash->{$items->[1]}->{reactions}->{$items->[0]} = 0;
	}
	$reactionhash->{$items->[0]}->{functions}->{$items->[1]} += $items->[2];
	$functionhash->{$items->[1]}->{reactions}->{$items->[0]} += $items->[2];
}
close($fhh);

open (my $out, ">", $directory."/ReactionRoleData.txt");
print $out "Reaction	Direction	Active	Essential	Active models	Essential models	Reactants with active side pathways	Reactants with essential side pathways	Average essential flux	Average active flux	Damage prone reactants	Role count	Max genomes	Max genome role	Max fusions	Max fusion role	Max fusion fraction	Max fusion fraction role	Max fusion score	Max fusion score role	Cpx	Transport	Models	DeltaG	DeltaG percentile\n";
foreach my $rxn (keys(%{$rxnacthash})) {
	print $out join("\t",@{$rxnacthash->{$rxn}})."\t";
	my $items = [split(/\t/,$rxn)];
	if (defined($reactionhash->{$items->[0]})) {
		$reactionhash->{$items->[0]}->{found} = 1;
		my $rolecount = keys(%{$reactionhash->{$items->[0]}->{functions}});
		my $maxgenomes = 0;
		my $maxgenomefunc;
		my $maxfusions = 0;
		my $maxfusionfunc;
		my $maxfusionfraction = 0;
		my $maxfusionfracfunc;
		my $maxfusionscore = 0;
		my $maxfusionscorefunc;
		foreach my $function (keys(%{$reactionhash->{$items->[0]}->{functions}})) {
			if ($function ne "hypotheticalprotein") {
				if (!defined($maxgenomefunc)) {
					$maxgenomes = $functionhash->{$function}->{genomes};
					$maxgenomefunc = $function;
					$maxfusions = $functionhash->{$function}->{fusiongenomes};
					$maxfusionfunc = $function;
					$maxfusionfraction = $functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes};
					$maxfusionfracfunc = $function;
					$maxfusionscore = $functionhash->{$function}->{fusiongenomes}*$functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes};
					$maxfusionscorefunc = $function;
				}
				if ($functionhash->{$function}->{genomes} > $maxgenomes) {
					$maxgenomefunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes} > $maxfusions) {
					$maxfusionfunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes} > $maxfusionfraction) {
					$maxfusionfracfunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes}*$functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes} > $maxfusionscore) {
					$maxfusionscorefunc = $function;
				}
			}
		}
		print $out $rolecount."\t".
			$maxgenomes."\t".
			$maxgenomefunc."\t".
			$maxfusions."\t".
			$maxfusionfunc."\t".
			$maxfusionfraction."\t".
			$maxfusionfracfunc."\t".
			$maxfusionscore."\t".
			$maxfusionscorefunc."\t".
			$reactionhash->{$items->[0]}->{cpx}."\t".
			$reactionhash->{$items->[0]}->{transport}."\t".
			$reactionhash->{$items->[0]}->{models};
	}
	print $out "\n";
}
foreach my $rxn (keys(%{$reactionhash})) {
	if (!defined($reactionhash->{$rxn}->{found})) {
		print $rxn."\n";
		print $out $rxn."											";
		$reactionhash->{$rxn}->{found} = 1;
		my $rolecount = keys(%{$reactionhash->{$rxn}->{functions}});
		my $maxgenomes = 0;
		my $maxgenomefunc;
		my $maxfusions = 0;
		my $maxfusionfunc;
		my $maxfusionfraction = 0;
		my $maxfusionfracfunc;
		my $maxfusionscore = 0;
		my $maxfusionscorefunc;
		foreach my $function (keys(%{$reactionhash->{$rxn}->{functions}})) {
			if ($function ne "hypotheticalprotein") {
				if (!defined($maxgenomefunc)) {
					$maxgenomes = $functionhash->{$function}->{genomes};
					$maxgenomefunc = $function;
					$maxfusions = $functionhash->{$function}->{fusiongenomes};
					$maxfusionfunc = $function;
					$maxfusionfraction = $functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes};
					$maxfusionfracfunc = $function;
					$maxfusionscore = $functionhash->{$function}->{fusiongenomes}*$functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes};
					$maxfusionscorefunc = $function;
				}
				if ($functionhash->{$function}->{genomes} > $maxgenomes) {
					$maxgenomefunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes} > $maxfusions) {
					$maxfusionfunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes} > $maxfusionfraction) {
					$maxfusionfracfunc = $function;
				}
				if ($functionhash->{$function}->{fusiongenomes}*$functionhash->{$function}->{fusiongenomes}/$functionhash->{$function}->{genomes} > $maxfusionscore) {
					$maxfusionscorefunc = $function;
				}
			}
		}
		print $out $rolecount."\t".
			$maxgenomes."\t".
			$maxgenomefunc."\t".
			$maxfusions."\t".
			$maxfusionfunc."\t".
			$maxfusionfraction."\t".
			$maxfusionfracfunc."\t".
			$maxfusionscore."\t".
			$maxfusionscorefunc."\t".
			$reactionhash->{$rxn}->{cpx}."\t".
			$reactionhash->{$rxn}->{transport}."\t".
			$reactionhash->{$rxn}->{models}."\n";
	}	
}
close($out);