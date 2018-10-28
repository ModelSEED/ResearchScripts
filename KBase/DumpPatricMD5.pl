#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use P3DataAPI;
use gjoseqlib;
use DateTime;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $path = $ARGV[0];

my $start = time();
my $failedcount = 0;
my $d = P3DataAPI->new();
my $genomes = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($path."/AllPatricGenomes.txt")}));
my $totalgenomes = keys(%{$genomes});
print "Total genomes:".$totalgenomes."\n";
my $failed = [];
for (my $i=0; $i < 100; $i++) { 
	open ( my $fh, ">", $path."/".$i."-genes.tsv") || Bio::KBase::ObjectAPI::utilities::error("Failure to open file, $!");
	print $fh "patric_id\tlocus_tag\tstart\tcontig\tproduct\tna_md5\taa_md5\n";
	my $count = 0;
	foreach my $genome (keys(%{$genomes})) {
		if ($genome =~ m/^(\d\d)/ || $genome =~ m/^(\d)\./) {
			my $storagekey = $1;
			if ($storagekey == $i) {
				#print "Processing ".$genome." (".$count.")\n";
				my(@fields) = qw(sequence_id start patric_id refseq_locus_tag product aa_sequence_md5 na_sequence_md5);
				my $ftrcount = 0;
				eval {
					my @res = $d->query("genome_feature",
						 ["select", @fields],
						 ["ne", "feature_type", "source"],
						 ["eq", "annotation", "PATRIC"],
						 ["eq", "genome_id", $genome],
					);
					print $fh "NEWGENOME:".$genome."\n";
					foreach my $ftr (@res) {
						if (defined($ftr->{na_sequence_md5})) {
							$ftrcount++;
							print $fh $ftr->{patric_id}."\t";
							if (defined($ftr->{refseq_locus_tag})) {
								print $fh $ftr->{refseq_locus_tag}
							}
							print $fh "\t".$ftr->{start}."\t".$ftr->{sequence_id}."\t";
							if (defined($ftr->{product})) {	
								print $fh $ftr->{product};
							}
							print $fh "\t".$ftr->{na_sequence_md5}."\t";
							if (defined($ftr->{aa_sequence_md5})) {
								print $fh $ftr->{aa_sequence_md5}
							}
							print $fh "\n";
						}
					}
				};
				if ($ftrcount == 0) {
					push(@{$failed},$genome);
					$failedcount++;
				} else {
					$count++;
				}
			}
		}
	}
	close ($fh);
	print $i."\t".$count."\n";
	print "Elapsed time:".(time()-$start)."\n";
	print "Failed count:".$failedcount."\n";
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/FailedGenomes.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($failed,1)]);