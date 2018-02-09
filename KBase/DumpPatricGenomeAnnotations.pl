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
my $output = {};
my $funchash = {};
my $funccount = 0;
my $failedcount = 0;
my $d = P3DataAPI->new();
my $genomes = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($path."/AllPatricGenomes.txt")}));
my $failed = [];
for (my $i=0; $i < 100; $i++) {
	my $output;
	my $count = 0;
	foreach my $genome (keys(%{$genomes})) {
		if ($genome =~ m/^(\d\d)/ || $genome =~ m/^(\d)\./) {
			my $storagekey = $1;
			if ($storagekey == $i) {
				#print "Processing ".$genome."\n";
				my(@fields) = qw(start product);
				eval {
					my @res = $d->query("genome_feature",
						 ["select", @fields],
						 ["ne", "feature_type", "source"],
						 ["eq", "annotation", "PATRIC"],
						 ["eq", "genome_id", $genome],
					);
					
					my $unsortedlist = [@res];
					foreach my $ftr (@{$unsortedlist}) {
						$ftr->{start} += 0;
					}
					my $sortedlist = [sort { $b->{start} cmp $a->{start} } @{$unsortedlist}];
					for my $ent (@{$sortedlist}) {
				    		next if !defined($ent->{product}) || $ent->{product} eq '';
					    	my $function = $ent->{product};
					  	my $array = [split(/\#/,$function)];
					  	$function = shift(@{$array});
						$function =~ s/\s+$//;
						$array = [split(/\s*;\s+|\s+[\@\/]\s+/,$function)];
						for (my $k=0; $k < @{$array}; $k++) {
							if (!defined($funchash->{$array->[$k]})) {
								$funchash->{$array->[$k]} = $funccount;
								$funccount++;
							}
							push(@{$output->{$genome}},$funchash->{$array->[$k]});
						}
					}
				};
				if (!defined($output->{$genome})) {
					push(@{$failed},$genome);
					$failedcount++;
				} else {
					$count++;
				}
			}
		}
	}
	if ($count > 0) {
		Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/".$i."-PatricAnnotations.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($output,1)]);
	}
	$output = {};
	print $i."\t".$count."\n";
	print "Function count:".$funccount."\n";
	print "Elapsed time:".(time()-$start)."\n";
	print "Failed count:".$failedcount."\n";
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/AllPatricFunctions.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($funchash,1)]);
Bio::KBase::ObjectAPI::utilities::PRINTFILE($path."/FailedGenomes.txt",[Bio::KBase::ObjectAPI::utilities::TOJSON($failed,1)]);