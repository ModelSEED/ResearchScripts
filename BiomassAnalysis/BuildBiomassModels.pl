use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;

my $procs = $ARGV[0];
my $index = $ARGV[1];

Bio::KBase::kbaseenv::create_context_from_client_config();
my $impl = fba_tools::fba_toolsImpl->new();

my $initial_genomes = Bio::KBase::kbaseenv::ws_client()->list_objects({
	workspaces => ["jplfaria:narrative_1492808527866"],
	type => "KBaseGenomes.Genome",
});
my $genomes;
my $genome_hash = {};
for (my $i=0; $i < @{$initial_genomes}; $i++) {
	if ($initial_genomes->[$i]->[1] =~ m/(.+)\.RAST$/) {
		$genome_hash->{$initial_genomes->[$i]->[1]} = $1;
		push(@{$genomes},$initial_genomes->[$i]->[1]);
	}
}

my $count = 0;
for (my $i=0; $i < @{$genomes}; $i++) {
	if ($i % $procs  == $index) {
		if ($count % 100 == 0) {
			$impl = fba_tools::fba_toolsImpl->new();		}
		eval {
			$impl->build_metabolic_model({
				workspace => "chenry:narrative_1493181437626",
				genome_id => $genomes->[$i],
				fbamodel_output_id => $genome_hash->{$genomes->[$i]}.".mdl",
				media_id => "Carbon-D-Glucose",
				template_id => "SuperBiomassTemplate",
				template_workspace => "chenry:narrative_1493181437626",
				genome_workspace => "jplfaria:narrative_1492808527866",
				media_workspace => "KBaseMedia",
				gapfill_model => 1
			});
		};
		if ($@) {
			my $error = $@;
			print $i.":".$genomes->[$i].":".$error;
		} else {
			print $i.":".$genomes->[$i].":success";
		}
		$count++;
	}
}