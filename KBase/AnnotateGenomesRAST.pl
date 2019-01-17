#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $ws = Bio::KBase::kbaseenv::ws_client();

my $workspace = $ARGV[0];
#$workspace = "filipeliu:narrative_1529802940504";
#$workspace = "chenry:narrative_1524167538737";

my $objects = Bio::KBase::kbaseenv::list_objects({
	workspaces => [$workspace],
	type => "KBaseGenomes.Genome"
});

my $client = Bio::KBase::kbaseenv::rast_client();
for (my $i=0; $i < @{$objects}; $i++) {
	if ($objects->[$i]->[1] !~ m/\.RAST/) {
		print $i.":".$objects->[$i]->[1]."\n";
		my $orig_genome = $handler->util_get_object($workspace."/".$objects->[$i]->[1],{raw => 1});
		#my $genome_obj = Bio::KBase::ObjectAPI::KBaseGenomes::Genome->new($orig_genome);
		#my $data = $genome_obj->serializeToDB();
		my $newgenome = $client->run_pipeline($orig_genome,{stages => Bio::KBase::constants::gene_annotation_pipeline()});
		my $genehash = {};
		for (my $j=0; $j < @{$newgenome->{features}}; $j++) {
			$genehash->{$newgenome->{features}->[$j]->{id}} = $newgenome->{features}->[$j];
		}
		for (my $k=0; $k < @{$orig_genome->{features}}; $k++) {
			if (defined($genehash->{$orig_genome->{features}->[$k]->{id}}->{function})) {
				$orig_genome->{features}->[$k]->{functions} = [split(/\s*;\s+|\s+[\@\/]\s+/,$genehash->{$orig_genome->{features}->[$k]->{id}}->{function})];
			}
		}
		delete $orig_genome->{genbank_handle_ref};
		$handler->util_save_object($orig_genome,$workspace."/".$objects->[$i]->[1].".RAST2",{hash => 1, type => "KBaseGenomes.Genome"});
	}
}