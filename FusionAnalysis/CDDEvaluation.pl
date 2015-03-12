#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $inputdir = $ARGV[1];

my $CDDData = {};
my $functions = {};
my $cddfunctions = {};
my $array;
open(my $fh, "<", $inputdir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	my $fh;
	open($fh, "<", $inputdir.$array->[$i]);
	my $line = <$fh>;
	my $correlate = {};
	while (my $line = <$fh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		my $maxalign = $items->[7]/$items->[1];
		if (!defined($CDDData->{$items->[8]})) {
			$CDDData->{$items->[8]} = {
				functions => {},
				maxalign => 0,
				genecount => 0,
				singlegenecount => 0,
				singlegenes => []
			};
		}
		if (!defined($functions->{$items->[6]}->{$items->[8]})) {
			$functions->{$items->[6]}->{$items->[8]} = 0;
			$cddfunctions->{$items->[8]}->{$items->[6]} = 0;
		}
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

open(GENELIST, ">", $directory."GeneList.txt");
open(CDDFUNCS, ">", $directory."CDDFunctions.txt");
open(CDDLIST, ">", $directory."CDDList.txt");
print GENELIST "CDD\tSingleGenes\n";
print CDDFUNCS "CDD\tFunctions\n";
print CDDLIST "CDD\tMax align\tSingleGenes\tTotalGenes\n";
foreach my $key (keys(%{$CDDData})) {
	print CDDFUNCS $key."\t";
	my $first = 0;
	foreach my $func (keys(%{$cddfunctions->{$key}})) {
		print CDDFUNCS $func.":::".$cddfunctions->{$key}->{$func};
		if ($first != 0) {
			print CDDFUNCS ";;;";
		}
		$first = 1;
	}
	print CDDFUNCS "\n";
	print GENELIST $key."\t".join(";",@{$CDDData->{$items->[8]}->{singlegenes}})."\n";
	print CDDLIST $key."\t".$CDDData->{$items->[8]}->{maxalign}."\t".$CDDData->{$items->[8]}->{singlegenecount}."\t".$CDDData->{$items->[8]}->{genecount}."\n";
}
close(CDDLIST);
close(GENELIST);
close(CDDFUNCS);

open(FUNCCDDS, ">", $directory."FunctionCDDs.txt");
print FUNCCDDS "Function\tCDDs\n";
foreach my $key (keys(%{$functions})) {
	print FUNCCDDS $key."\t";
	my $first = 0;
	foreach my $cdd (keys(%{$functions->{$key}})) {
		print FUNCCDDS $cdd.":::".$functions->{$key}->{$cdd};
		if ($first != 0) {
			print FUNCCDDS ";;;";
		}
		$first = 1;
	}
	print FUNCCDDS "\n";
}
close(FUNCCDDS);

1;