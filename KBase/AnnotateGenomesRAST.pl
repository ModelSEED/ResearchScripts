#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $ws = Bio::KBase::kbaseenv::ws_client();

my $workspace = "jplfaria:narrative_1524466549180";
#my $workspace = "chenry:narrative_1524167538737";

my $objects = Bio::KBase::kbaseenv::list_objects({
	workspaces => [$workspace],
	type => "KBaseGenomes.Genome"
});

my $client = Bio::KBase::kbaseenv::rast_client();
for (my $i=0; $i < @{$objects}; $i++) {
	if ($objects->[$i]->[1] !~ m/\.RAST/) {
		print $i.":".$objects->[$i]->[1]."\n";
		my $genome = $handler->util_get_object($workspace."/".$objects->[$i]->[1]);
		my $data = $genome->serializeToDB();
		my $newgenome = $client->run_pipeline($data,{stages => Bio::KBase::constants::gene_annotation_pipeline()});
		my $genehash = {};
		for (my $j=0; $j < @{$newgenome->{features}}; $j++) {
			#print $newgenome->{features}->[$j]->{id}."\t".$newgenome->{features}->[$j]->{function}."\n";
			$genehash->{$newgenome->{features}->[$j]->{id}} = $newgenome->{features}->[$j];
		}
		my $ftrs = $genome->features();
		for (my $k=0; $k < @{$ftrs}; $k++) {
			if (defined($genehash->{$ftrs->[$k]->id()}->{function})) {
				$ftrs->[$k]->function($genehash->{$ftrs->[$k]->id()}->{function});
			}
		}
		$handler->util_save_object($genome,$workspace."/".$objects->[$i]->[1].".RAST2");
	}
}