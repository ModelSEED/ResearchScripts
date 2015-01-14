#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];

my $objects = get_ws_objects_list($workspace);
my $hash;
for (my $i=0;$i < @{$objects}; $i++) {
	if (!defined($hash->{$objects->[$i]->[2]})) {
		$hash->{$objects->[$i]->[2]} = 0;
	}
	$hash->{$objects->[$i]->[2]}++;
}

foreach my $type (keys(%{$hash})) {
	print $type."\t".$hash->{$type}."\n";
}