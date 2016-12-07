#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

my $ws = "chenry:1431835409789";
(my $data,my $prov) = get_workspace_object("chenry:1432528044619/Homo_sapien_Genome");
my $ftrs = $data->{features};
my $ftrhash;
my $origid = [];
for (my $i=0; $i < @{$ftrs}; $i++) {
	$origid->[$i] = $ftrs->[$i]->{id};
	$ftrhash->{$ftrs->[$i]->{location}->[0]->[0].":".$ftrs->[$i]->{id}} = $ftrs->[$i];
}
my ($contig,$protein,$gene);
open(my $fh, "<", $directory."/genome.gbk");
my $done;
my $notfound = 0;
my $repeat = 0;
my $missed = 0;
while (my $line = <$fh>) {
	chomp($line);
	if ($line =~ m/LOCUS\s+(\w+)\s+/) {
		if (defined($contig) && defined($protein) && defined($gene)) {
			if (!defined($ftrhash->{$contig.":".$protein})) {
				$notfound++;
				print $contig.":".$protein." not found!\n";
			} else {
				if (defined($done->{$contig.":".$protein})) {
					$repeat++;
					print "Repeat ".$done->{$contig.":".$protein}."\t".$gene."\n";
				}
				$done->{$contig.":".$protein} = $gene;
				$ftrhash->{$contig.":".$protein}->{id} = $gene;
				$ftrhash->{$contig.":".$protein}->{aliases} = [$contig.":".$protein];
			}
		}
		$contig = $1;
	} elsif ($line =~ m/\/gene=\"(\w+)\"/) {
		$protein = $1;
	} elsif ($line =~ m/\/db_xref=\"CCDS:CCDS([\w\.]+)\"/) {
		$gene = $1;
	}
}
close($fh);

my $good = 0;
for (my $i=0; $i < @{$ftrs}; $i++) {
	if (!defined($done->{$ftrs->[$i]->{location}->[0]->[0].":".$origid->[$i]})) {
		$missed++;
		print "Missed:".$ftrs->[$i]->{location}->[0]->[0].":".$ftrs->[$i]->{id}."\n";
	} else {
		$good++;
	}
}
print "Not found:".$notfound."\n";
print "Repeat:".$repeat."\n";
print "Missed:".$missed."\n";
print "Good:".$good."\n";

$data->{contigset_ref} = "PublishedModelGenomes/Homo_sapien.contigs";
save_workspace_object("PublishedModelGenomes/Homo_sapien",$data,"KBaseGenomes.Genome");
