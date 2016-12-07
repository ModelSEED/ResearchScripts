#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $genome = $ARGV[0];
my $filename = $ARGV[1];

(my $obj,my $info) = get_workspace_object($genome);

open INPUT, ">", $filename or do {
	die "$0: open ".$filename.": $!";
};
my $ftrs = $obj->{features};
for (my $i=0; $i < @{$ftrs}; $i++) {
	if ($ftrs->[$i]->{protein_translation}) {
		print INPUT ">".$ftrs->[$i]->{id}."\n".$ftrs->[$i]->{protein_translation}."\n";
	}
}
close(INPUT);

