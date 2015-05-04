#!/usr/bin/perl -w

use strict;
use Config::Simple;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
use File::Path;
$|=1;

my $directory = $ARGV[0];
my $procindex = $ARGV[1];
my $numprocs = $ARGV[2];

my $list = get_ws_objects_list("KBasePublicModelsV4","KBaseFBA.FBAModel");
my $mdllist;
for (my $i=0; $i < @{$list}; $i++) {
	if ($list->[$i]->[1] =~ m/\.gf$/) {
		push(@{$mdllist},$list->[$i]);
	}
}

my $fba = get_fba_client();
for (my $i=0; $i < @{$mdllist}; $i++) {
	my $value = $i-$procindex;
	if (($value % $numprocs) == 0) {
		my $sbml = $fba->export_fbamodel({
			model => $mdllist->[$i]->[1],
			workspace => "KBasePublicModelsV4",
			format => "sbml"
		});
		if ($sbml =~ m/model\sid.+name=\"(.+)_SEED_model/) {
			my $name = $1;
			$name =~ s/\s/_/g;
			open ( my $fh, ">", $directory."/".$name.".xml");
			print $fh $sbml;
			close($fh);
		}
	}
}