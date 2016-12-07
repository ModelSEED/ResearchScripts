#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use FindBin qw($Bin);
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $ws = "chenry:1431835409789";

(my $data,my $prov) = get_workspace_object($ws."/Bacillus_subtilis_DSM_SEED");
my $ftrhash = {};
for (my $i=0; $i < @{$data->{features}}; $i++) {
	$ftrhash->{$data->{features}->[$i]->{id}} = $data->{features}->[$i];
}

my $maxvalues = {};
my $expdata = {};
print $Bin."/TranscriptomicsRaw.txt"."\n";
open(my $fh, "<", $Bin."/TranscriptomicsRaw.txt");
my $line = <$fh>;
chomp($line);
my $headings = [split(/\t/,$line)];
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if ($array->[0] =~ m/BSU(.+)/) {
		my $newid = "BSU6051_".$1;
		if (defined($ftrhash->{$newid})) {
			for (my $j=1; $j < @{$headings}; $j++) {
				if (!defined($maxvalues->{$headings->[$j]})) {
					$maxvalues->{$headings->[$j]} = 0;
				}
				if ($array->[$j] > $maxvalues->{$headings->[$j]}) {
					$maxvalues->{$headings->[$j]} = $array->[$j];	
				}
				$expdata->{$headings->[$j]}->{$newid} = $array->[$j];
			}
		} else {
			print "Could not find ".$newid."\n";
		}
	}
}
close($fh);

foreach my $run (keys(%{$expdata})) {
	print $run."\t".$maxvalues->{$run}."\n";
	open(my $fhh, ">", $Bin."/".$run.".txt");
	foreach my $gene (keys(%{$expdata->{$run}})) {
		print $fhh $gene."\t".$expdata->{$run}->{$gene}/$maxvalues->{$run}."\n";
	}
	close($fhh);
}