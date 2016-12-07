#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

(my $obj,my $info) = get_workspace_object("10559/Escherichia_coli_K_12_substr_MG1655.kbase");

my $ftrs = $obj->{features};
print "id\ttype\tfunction\tstart\tlength\tdirection\tprotein\n";
for (my $i=0; $i < @{$ftrs}; $i++) {
	print $ftrs->[$i]->{id}."\t".
		$ftrs->[$i]->{type}."\t".
		$ftrs->[$i]->{function}."\t".
		$ftrs->[$i]->{location}->[0]->[0]."\t".
		$ftrs->[$i]->{location}->[0]->[1]."\t".
		$ftrs->[$i]->{location}->[0]->[3]."\t".
		$ftrs->[$i]->{location}->[0]->[2]."\t".
		$ftrs->[$i]->{protein_translation}."\n";
}