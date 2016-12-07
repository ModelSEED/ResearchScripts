#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

(my $object,my $info) = get_workspace_object("chenry:1435635808596/Klebsiella_pneumoniae_kppr1","KBaseFBA.FBA");

my $ftrs = $object->{features};
my $contigftrs = {};
foreach my $ftr (@{$ftrs}) {
	push(@{$contigftrs->{$ftr->{location}->[0]->[0]}},$ftr);
}
my $contigoperons = {};
foreach my $contig (keys(%{$contigftrs})) {
	my $contigftrs->{$contig} = [sort { $a->{location}->[0]->[1] >= $b->{location}->[0]->[1] } @{$contigftrs->{$contig}}];
	my $operon;
	my $currdir;
	my $laststop;
	foreach my $ftr (@{$contigftrs->{$contig}}) {
		my $truestart = $ftr->{location}->[0]->[1];
		my $truestop = $truestart+$ftr->{location}->[0]->[3];
		if ($ftr->{location}->[0]->[2] eq "-") {
			$truestop = $truestart;
			$truestart = $truestart-$ftr->{location}->[0]->[3];
		}
		if (defined($currdir)) {
			if ($currdir ne $ftr->{location}->[0]->[2] || ($truestart-$laststop) > 200) {
				push(@{$contigoperons->{$contig}},$operon);
				$operon = [$ftr];
			} elsif ($currdir eq "-") {
				unshift(@{$operon},$ftr);
			} elsif ($currdir eq "+") {
				push(@{$operon},$ftr);
			}
		}
		$currdir = $ftr->{location}->[0]->[2];
		$laststop = $truestop;
	}
	push(@{$contigoperons->{$contig}},$operon);
}
my $kotranslation;

#print $out $ftr->{id}."\t".$ftr->{location}->[0]->[0]."\t".$ftr->{location}->[0]->[1]."\t".$ftr->{location}->[0]->[2]."\t".$ftr->{location}->[0]->[3]."\t".$ftr->{function}."\n";
#open (my $out, ">", $directory."/KlebsiellaGenome.txt");
#print $out "ID\tContig\tStart\tDirection\tLength\tFunction\n";
#close($out);