#!/usr/bin/perl -w

use strict;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
$|=1;

my $directory = $ARGV[0];
print $directory."KBaseGenomeList.txt"."\n";
open ( MAIN, ">", $directory."KBaseGenomeList.txt");
my $list = get_ws_objects_list("KBasePublicGenomesV5","KBaseGenomes.Genome");
print Data::Dumper->Dump([$list->[100]]);
my $headers = [keys(%{$list->[100]->[10]})];
print MAIN "id";
foreach my $key (keys(%{$list->[100]->[10]})) {
	print MAIN "\t".$key;
}
print MAIN "\n";
for (my $i=0; $i < @{$list}; $i++) {
	print MAIN $list->[$i]->[1];
	foreach my $key (@{$headers}) {
		print MAIN "\t";
		if (defined($list->[$i]->[10]->{$key})) {
			print MAIN $list->[$i]->[10]->{$key};
		}
	}
	print MAIN "\n";
}
close(MAIN);