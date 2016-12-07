#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $model = $ARGV[0];

(my $obj,my $info) = get_workspace_object($model);

my $rxns = $obj->{modelreactions};
my $genehash;
for (my $i=0; $i < @{$rxns}; $i++) {
	my $prots = $rxns->[$i]->{modelReactionProteins};
	my $exclusive = 1;
	if (@{$prots} > 1) {
		$exclusive = 0;
	}
	for (my $j=0; $j < @{$prots}; $j++) {
		my $subunits = 	$prots->[$j]->{modelReactionProteinSubunits};
		for (my $k=0; $k < @{$subunits}; $k++) {
			my $ftrs = $subunits->[$k]->{feature_refs};
			if (@{$ftrs} > 1) {
				$exclusive = 0;
			}
			for (my $m=0; $m < @{$ftrs}; $m++) {
				if ($ftrs->[$m] =~ m/\/([^\/]+)$/) {
					push(@{$genehash->{$1}},$rxns->[$i]->{id}.":".$exclusive);
				}	
			}
		}
	}
}

print "id\treactions\n";
foreach my $gene (keys(%{$genehash})) {
	print $gene."\t".join("|",@{$genehash->{$gene}})."\n";
}