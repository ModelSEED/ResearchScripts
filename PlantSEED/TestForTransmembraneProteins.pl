#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;
use Data::Dumper;
use Getopt::Long;
use ScriptThing;
use LWP;
use URI;

$| = 1;

my $directory = $ARGV[0];

my $sapsvr = ModelSEED::Client::SAP->new();

print "Loading transporter candidates!\n";
open(my $fhhh, "<", $directory."Candidates.txt");
my $genehash;
while (my $line = <$fhhh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$genehash->{$array->[0]} = [split(/;/,$array->[1])];
}
close($fhhh);

my $genelist = [keys(%{$genehash})];
my $geneseqs = $sapsvr->ids_to_sequences(
	-ids => $genelist,
	-protein => 1,
	-fasta   => 0
);

print "Gene\tInside helices\tOutside helices\tModel\tCount\n";
foreach my $gene (keys(%{$geneseqs})) {
	my $seq = $geneseqs->{$gene};
	my $browser = LWP::UserAgent->new;
	my $url = "http://www.ch.embnet.org/cgi-bin/TMPRED_form_parser";
	my $response = $browser->post($url,[
		'outmode' => "html",
		'min' => "17",
		'max' => "33",
		'comm' => "FIGSEQ",
		'format' => "plain_text",
		'seq' => $seq
	]);
	my $x = $response->content;
	my $array = [split(/\n/,$x)];
	my $insidecount = 0;
	my $outsidecount = 0;
	my $model = 1;
	my $count = @{$genehash->{$gene}};
	for (my $i=0; $i < @{$array}; $i++) {
		if ($array->[$i] =~ m/Inside\sto\soutside\shelices\s:\s+(\d+)\sfound/) {
			my $count = $1;
			$i++;
			for (my $j=0; $j < $count; $j++) {
				$i++;
				if ($array->[$i] =~ m/\s*(\d+)[\(\s]+(\d+)[\)\s]+(\d+)[\(\s]+(\d+)[\)\s]+(\d+)\s+(\d+)/) {
					if ($5 >= 500) {
						$insidecount++;
					}
				}
			}
		}
		if ($array->[$i] =~ m/Outside\sto\sinside\shelices\s:\s+(\d+)\sfound/) {
			my $count = $1;
			$i++;
			for (my $j=0; $j < $count; $j++) {
				$i++;
				if ($array->[$i] =~ m/\s*(\d+)[\(\s]+(\d+)[\)\s]+(\d+)[\(\s]+(\d+)[\)\s]+(\d+)\s+(\d+)/) {
					if ($5 > 500) {
						$outsidecount++;
					}
				}
			}
		}
		if ($array->[$i] =~ m/probably\sno\stransmembrane\sprotein/) {
			$model = 0;
		}
	}
	print $gene."\t".$insidecount."\t".$outsidecount."\t".$model."\t".$count."\n";
}