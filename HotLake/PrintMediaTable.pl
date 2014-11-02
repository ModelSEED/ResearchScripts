#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];

(my $obj,my $info) = get_workspace_object("KBasePhenotypeDatasets/83333.1-biolog");

my $CompoundHash;
my $MediaHash = {};
my $phenos = $obj->{phenotypes};
for (my $i=0; $i < @{$phenos}; $i++) {
	my $pheno = $phenos->[$i];
	(my $media,$info) = get_workspace_object($pheno->{media_ref});
	my $mediacpds = $media->{mediacompounds};
	for (my $j=0; $j < @{$mediacpds}; $j++) {
		my $cpd = $mediacpds->[$j]->{compound_ref};
		if ($cpd =~ m/(cpd\d+)/) {
			$CompoundHash->{$1} = "";
			$MediaHash->{$media->{name}}->{$1} = 1;
		}
	}
}

my $fba = get_fba_client();
my $output = $fba->get_compounds({
	compounds => [keys(%{$CompoundHash})]
});
for (my $i=0; $i < @{$output}; $i++) {
	$CompoundHash->{$output->[$i]->{id}} = $output->[$i]->{name};
}

foreach my $media (keys(%{$MediaHash})) {
	print $media."\t";
	foreach my $id (keys(%{$MediaHash->{$media}})) {
		print $CompoundHash->{$id}.";"
	}
	print "\n";
}