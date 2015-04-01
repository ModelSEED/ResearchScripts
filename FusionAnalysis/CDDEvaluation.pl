#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];
my $genome = $ARGV[2];
my $CDDData = {};
my $cddfunctions = {};
my $array;
open(my $fh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);
if (defined($genome)) {
	$array = [$genome];
}
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $correlate = {};
	my $genes = {};
	while (my $line = <$fh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		my $maxalign = 3*$items->[7]/$items->[1];
		if (!defined($CDDData->{$items->[8]})) {
			$CDDData->{$items->[8]} = {
				functions => {},
				maxalign => 0,
				genecount => 0,
				singlegenecount => 0,
				singlegenes => [],
				size => 0
			};
		}
		if ($items->[11] > $CDDData->{$items->[8]}->{size}) {
			$CDDData->{$items->[8]}->{size} = $items->[11];
		}
		if (!defined($functions->{$items->[6]}->{$items->[8]})) {
			$functions->{$items->[6]}->{$items->[8]} = 0;
			$cddfunctions->{$items->[8]}->{$items->[6]} = 0;
		}
		$genes->{$items->[0]}->{$items->[8]} = [$items->[3],$items->[4]];
		$functions->{$items->[6]}->{$items->[8]} += $maxalign*$items->[5];
		$cddfunctions->{$items->[8]}->{$items->[6]} += $maxalign*$items->[5];
		$CDDData->{$items->[8]}->{genecount}++;
		if ($maxalign > 0.9) {
			$CDDData->{$items->[8]}->{singlegenecount}++;
			push(@{$CDDData->{$items->[8]}->{singlegenes}},$items->[13]);
		}
		if ($maxalign > $CDDData->{$items->[8]}->{maxalign}) {
			$CDDData->{$items->[8]}->{maxalign} = $maxalign;
		}
	}
}

my $functions = {};
open(GENELIST, ">", $directory."GeneList.txt");
open(CDDFUNCS, ">", $directory."CDDFunctions.txt");
open(CDDLIST, ">", $directory."CDDList.txt");
print GENELIST "CDD\tSingleGenes\n";
print CDDFUNCS "CDD\tFunctions\n";
print CDDLIST "CDD\tMax align\tSize\tSingleGenes\tTotalGenes\n";
foreach my $key (keys(%{$CDDData})) {
	print CDDFUNCS $key."\t";
	my $first = 0;
	my $array = [];
	my $sum = 0;
	foreach my $func (keys(%{$cddfunctions->{$key}})) {
		$sum += $cddfunctions->{$key}->{$func};
		push(@{$array},[$func,$cddfunctions->{$key}->{$func}]);
	}
	for (my $i=0; $i < @{$array}; $i++) {
		$array->[$i]->[1] = $array->[$i]->[1]/$sum;
		$cddfunctions->{$key}->{$array->[$i]->[0]} = $array->[$i]->[1];
		$functions->{$array->[$i]->[0]}->{$key} = $array->[$i]->[1];
	}
	my $sortedarray = [sort { $b->[1] <=> $a->[1] } @{$array}];
	for (my $i=0; $i < @{$sortedarray}; $i++) {
		if ($first != 0) {
			print CDDFUNCS ";;;";
		}
		print CDDFUNCS $sortedarray->[$i]->[0].":::".$sortedarray->[$i]->[1];
		$first = 1;
	}
	print CDDFUNCS "\n";
	print GENELIST $key."\t".join(";",@{$CDDData->{$key}->{singlegenes}})."\n";
	print CDDLIST $key."\t".$CDDData->{$key}->{maxalign}."\t".$CDDData->{$key}->{size}."\t".$CDDData->{$key}->{singlegenecount}."\t".$CDDData->{$key}->{genecount}."\n";
}
close(CDDLIST);
close(GENELIST);
close(CDDFUNCS);

open(FUNCCDDS, ">", $directory."FunctionCDDs.txt");
print FUNCCDDS "Function\tCDDs\n";
foreach my $key (keys(%{$functions})) {
	print FUNCCDDS $key."\t";
	my $first = 0;
	my $array = [];
	my $sum = 0;
	foreach my $cdd (keys(%{$functions->{$key}})) {
		$sum += $functions->{$key}->{$cdd};
		push(@{$array},[$cdd,$functions->{$key}->{$cdd}]);
	}
	for (my $i=0; $i < @{$array}; $i++) {
		$array->[$i]->[1] = $array->[$i]->[1]/$sum;
		$functions->{$key}->{$array->[$i]->[0]} = $array->[$i]->[1];
	}
	my $sortedarray = [sort { $b->[1] <=> $a->[1] } @{$array}];
	for (my $i=0; $i < @{$sortedarray}; $i++) {
		if ($first != 0) {
			print FUNCCDDS ";;;";
		}
		print FUNCCDDS $sortedarray->[$i]->[0].":::".$sortedarray->[$i]->[1];
		$first = 1;
	}
	print FUNCCDDS "\n";
}
close(FUNCCDDS);

1;