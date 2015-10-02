#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

(my $obj,my $info) = get_workspace_object("chenry:1435635808596/Klebsiella_pneumoniae_kppr1");
(my $objseed,my $info) = get_workspace_object("chenry:1435635808596/Klebsiella_pneumoniae_kppr1_Reanno");

my $ftrs = $obj->{features};
print "id\ttype\tfunction\tSEED function\tstart\tlength\tdirection\tprotein\n";
for (my $i=0; $i < @{$ftrs}; $i++) {
	print $ftrs->[$i]->{id}."\t".
		$ftrs->[$i]->{type}."\t".
		$ftrs->[$i]->{function}."\t".
		$objseed->{features}->[$i]->{function}."\t".
		$ftrs->[$i]->{location}->[0]->[0]."\t".
		$ftrs->[$i]->{location}->[0]->[1]."\t".
		$ftrs->[$i]->{location}->[0]->[3]."\t".
		$ftrs->[$i]->{location}->[0]->[2]."\t".
		$ftrs->[$i]->{protein_translation}."\n";
}