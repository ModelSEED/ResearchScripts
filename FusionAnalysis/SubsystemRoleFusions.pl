#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];

my $subsysclasshash;
my $funchash = {};
open(my $fh, "<", $directory."FunctionFusionAnalysis.txt");
my $line = <$fh>;
while ($line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$funchash->{$array->[0]} = {
		searchname => $array->[0],
		count => $array->[1],
		genomecount => $array->[2],
		fusions => $array->[3],
		fusiongenomes => $array->[4],
		fractionfusions => 0,
		fractionfusiongenomes => 0,
		score => 0,
		subsystems => {},
		subsystem => "none",
		classone => "none",
		classtwo => "none",
		subsysfusfunc => 0,
		subsysfusfuncgen => 0
	};
	if ($array->[1] > 0) {
		$funchash->{$array->[0]}->{fractionfusions} = $array->[3]/$array->[1];
	}
	if ($array->[2] > 0) {
		$funchash->{$array->[0]}->{fractionfusiongenomes} = $array->[4]/$array->[2];
		$funchash->{$array->[0]}->{score} = $array->[4]*$array->[4]/$array->[2];
	}
}
close($fh);
my $sapsvr = ModelSEED::Client::SAP->new();
my $data = $sapsvr->all_subsystems({
	-usable => 1,
	-exclude => ["cluster-based", "experimental"],
	-aux => 1
});
my $sshash = {};
foreach my $subsys (keys(%{$data})) {
	my $funccount = @{$data->{$subsys}->[2]};
	$sshash->{$subsys} = {
		name => $subsys,
		class => $data->{$subsys}->[1]->[0],
		subclass => $data->{$subsys}->[1]->[1],
		functions => $funccount,
		fusionfunc => 0,
		fusiongenomes => 0,
		genomes => 0,
		functionhash => {}
	};
	if (!defined($subsysclasshash->{$data->{$subsys}->[1]->[0]})) {
		$subsysclasshash->{$data->{$subsys}->[1]->[0]} = {
			subsystems => 0,
			functions => 0,
			fusionfunc => 0,
			fusiongenomes => 0,
			genomes => 0,
			functionhash => {}
		};
	}
	$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{subsystems}++;
	foreach my $func (@{$data->{$subsys}->[2]}) {
		$func = lc($func);
		$func =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
		$func =~ s/\s//g;
		$func =~ s/\#.*$//g;
		$func =~ s/\(ec\)//g;
		if (!defined($subsysclasshash->{$data->{$subsys}->[1]->[0]}->{functionhash}->{$func})) {
			$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{functionhash}->{$func} = 1;
			$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{genomes} += $funchash->{$func}->{genomecount};
			$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{fusiongenomes} += $funchash->{$func}->{fusiongenomes};
			$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{functions}++;
			if ($funchash->{$func}->{fractionfusiongenomes} >= 0.611) {
				$subsysclasshash->{$data->{$subsys}->[1]->[0]}->{fusionfunc}++
			}
		}
		$sshash->{$subsys}->{functionhash}->{$func} = 1;
		$sshash->{$subsys}->{genomes} += $funchash->{$func}->{genomecount};
		$sshash->{$subsys}->{fusiongenomes} += $funchash->{$func}->{fusiongenomes};
		if ($funchash->{$func}->{fractionfusiongenomes} >= 0.611) {
			$sshash->{$subsys}->{fusionfunc}++;
		}
		$funchash->{$func}->{subsystems}->{$subsys} = 1;
	}
}

foreach my $key (keys(%{$funchash})) {
	my $maxfusion = 0;
	my $subsys;
	foreach my $ss (keys(%{$funchash->{$key}->{subsystems}})) {
		if (!defined($subsys) || $sshash->{$ss}->{fusiongenomes} > $maxfusion) {
			$maxfusion = $sshash->{$ss}->{fusiongenomes};
			$subsys = $ss;
		}
	}
	if (defined($subsys)) {
		$funchash->{$key}->{subsystem} = $sshash->{$subsys}->{name};
		$funchash->{$key}->{classone} = $sshash->{$subsys}->{class};
		$funchash->{$key}->{classtwo} = $sshash->{$subsys}->{subclass};
		$funchash->{$key}->{subsysfusfunc} = $sshash->{$subsys}->{fusionfunc};
		$funchash->{$key}->{subsysfusfuncgen} = $sshash->{$subsys}->{fusiongenomes};
	}
}
open (my $out, ">", $directory."/FunctionFusionAnalysisTwo.txt");
print $out "Function\tCount\tGenome count\tFusions\tFusion genomes\tFraction fusions\tFraction fusion genomes\tScore\tSubsystem\tClass one\tClass two\tFusion functions\tFunction fusion genomes\n";
foreach my $key (keys(%{$funchash})) {
	if ($funchash->{$key}->{genomecount} >= 10 && length($funchash->{$key}->{searchname}) > 0 && $funchash->{$key}->{searchname} !~ m/predicted|hypothetical|putative|possible|probable/) {
		print $out $funchash->{$key}->{searchname}."\t".
			$funchash->{$key}->{count}."\t".
			$funchash->{$key}->{genomecount}."\t".
			$funchash->{$key}->{fusions}."\t".
			$funchash->{$key}->{fusiongenomes}."\t".
			$funchash->{$key}->{fractionfusions}."\t".
			$funchash->{$key}->{fractionfusiongenomes}."\t".
			$funchash->{$key}->{score}."\t".
			$funchash->{$key}->{subsystem}."\t".
			$funchash->{$key}->{classone}."\t".
			$funchash->{$key}->{classtwo}."\t".
			$funchash->{$key}->{subsysfusfunc}."\t".
			$funchash->{$key}->{subsysfusfuncgen}."\n";
	}
}
close($out);
open (my $outt, ">", $directory."/SubsystemFusions.txt");
print $outt "Subsystem\tClass one\tClass two\tType\tFunctions\tFusion functions\tFraction fusions\tFunction fusion genomes\tFunction genomes\tFraction fusion genomes\n";
foreach my $key (keys(%{$sshash})) {
	if ($sshash->{$key}->{functions} > 0 && $sshash->{$key}->{genomes} > 0) {
		my $funcfraction = $sshash->{$key}->{fusionfunc}/$sshash->{$key}->{functions};
		my $funcgenfraction = $sshash->{$key}->{fusiongenomes}/$sshash->{$key}->{genomes};
		print $outt $sshash->{$key}->{name}."\t".$sshash->{$key}->{class}."\t".$sshash->{$key}->{subclass}."\t".$sshash->{$key}->{type}."\t".$sshash->{$key}->{functions}."\t".$sshash->{$key}->{fusionfunc}."\t".$funcfraction."\t".$sshash->{$key}->{fusiongenomes}."\t".$sshash->{$key}->{genomes}."\t".$funcgenfraction."\n";
	}
}
close($outt);
open (my $outtt, ">", $directory."/SubsystemClassFusions.txt");
print $outtt "Class\tSubsystems\tFunctions\tFusion functions\tFraction fusions functions\tGenomes\tFusion genomes\tFraction fusion genomes\n";
foreach my $key (keys(%{$subsysclasshash})) {
	my $funcfraction = $subsysclasshash->{$key}->{fusionfunc}/$subsysclasshash->{$key}->{functions};
	my $funcgenfraction = $subsysclasshash->{$key}->{fusiongenomes}/$subsysclasshash->{$key}->{genomes};
	print $outtt $key."\t".$subsysclasshash->{$key}->{subsystems}."\t".$subsysclasshash->{$key}->{functions}."\t".$subsysclasshash->{$key}->{fusionfunc}."\t".$funcfraction."\t".$subsysclasshash->{$key}->{genomes}."\t".$subsysclasshash->{$key}->{fusiongenomes}."\t".$funcgenfraction."\n";
}
close($outtt);