#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];
my $genome = $ARGV[1];
my $workspace = $ARGV[2];

(my $obj,my $info) = get_workspace_object($workspace."/".$genome);
my $proteins;
my $output = "";
my $ftrs = $obj->{features};
for (my $i=0; $i < @{$ftrs}; $i++) {
	if ($ftrs->[$i]->{protein_translation}) {
		$output .= $ftrs->[$i]->{id}."\t".3*length($ftrs->[$i]->{protein_translation})."\n";
	}
}
open(my $fh, ">", $directory.$genome.".fasta");
print $fh $output;
close ($fh);