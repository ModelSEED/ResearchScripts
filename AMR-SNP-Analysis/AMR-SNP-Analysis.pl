#!/usr/bin/perl
$|=1;
my $directory = $ARGV[0];

my $data;
open(my $fh, "<", $directory."rawdata.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$data},[split(/\t/,$line)]);
}
close($fh);

open($fh, ">", $directory."amr_snp.arff");
print $fh "\@relation amr-snps\n\n";
print $fh "\@attribute Strain string\n";
print $fh "\@attribute Phenotype {1,2,3,4}\n";
for (my $i=2; $i < @{$data->[0]}; $i++) {
	my $values = {};
	for (my $j=0; $j < @{$data}; $j++) {
		$values->{$data->[$j]->[$i]} = 1;
	}
	print $fh "\@attribute SNP".($i-1)." {".join(",",keys(%{$values}))."}\n";
}
print $fh "\n\@data\n";
for (my $i=0; $i < @{$data}; $i++) {
	print $fh join(",",@{$data->[$i]})."\n";
}
close($fh);

1;