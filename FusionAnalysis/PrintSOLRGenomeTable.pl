#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use ModelSEED::Client::SAP;

my $directory = $ARGV[0];

my $genomehash;
open(my $fh, "<", $directory."/TableS9-GenomeStats.txt");
my $line = <$fh>;
while ($line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	$genomehash->{$items->[0]} = $line;
}
close($fh);

my $list = [keys(%{$genomehash})];
for (my $i=0; $i < @{$list}; $i++) {
	$list->[$i] =~ s/fig\|//;
}

my $sapsvr = ModelSEED::Client::SAP->new();
my $headings = ["dna-size","contigs","complete","domain","gc-content","taxonomy","md5","rnas","genetic-code"];
my $genomedata = $sapsvr->genome_data({
	-ids => $list,
	-data => $headings
});

for (my $i=0; $i < @{$list}; $i++) {
	$genomedata->{"fig|".$list->[$i]} = $genomedata->{$list->[$i]};
}

open (my $out, ">", $directory."/../SOLR/SOLR-GenomeStats.txt");
print $out "id\tname\tgenes\tgenes_with_cdds\ttotal_cdd_hits\tgenes_with_nonoverlapping_cdds\tfull_genes_cdd_hits\tgenes_with_nonoverlapping_full_gene_cdds\tfinal_predicted_fusions\tfraction_final_fusion_predictions\t".join("\t",@{$headings})."\n";
foreach my $key (keys(%{$genomehash})) {
	print $out $genomehash->{$key}."\t";
	if (defined($genomedata->{$key})) {
		print $out join("\t",@{$genomedata->{$key}})."\n";
	} else {
		print $out "\t\t\t\t\t\t\t\t\n";
	}
}
close($out);