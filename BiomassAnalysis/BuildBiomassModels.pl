use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;

my $procs = $ARGV[0];
my $index = $ARGV[1];

Bio::KBase::kbaseenv::create_context_from_client_config();
my $impl = fba_tools::fba_toolsImpl->new();

my $models = Bio::KBase::kbaseenv::ws_client()->list_objects({
	workspaces => ["jplfaria:narrative_1493091104031"],
	type => "KBaseFBA.FBAModel",
});
my $gramhash = {};
for (my $i=0; $i < @{$models}; $i++) {
	if ($models->[$i]->[1] =~ m/(.+)\.RAST\.mdl/) {
		my $genomeid = $1;
		my $model = $impl->util_get_object("jplfaria:narrative_1493091104031/".$models->[$i]->[1]);
		if ($model->biomasses()->[0]->name() =~ m/GramPositiveBiomass/) {
			$gramhash->{$genomeid} = "p";
		} else {
			$gramhash->{$genomeid} = "n";
		}
	}
}

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
			$impl = fba_tools::fba_toolsImpl->new();
		}
		eval {
			my $genomeid = $genome_hash->{$genomes->[$i]};
			my $tid;
			if (defined($gramhash->{$genomeid})) {
				if ($gramhash->{$genomeid} eq "p") {
					$tid = "GramPosModelTemplate";
				} else {
					$tid = "GramNegModelTemplate";
				}
				$impl->build_metabolic_model({
					workspace => "chenry:narrative_1493799658132",
					genome_id => $genomes->[$i],
					fbamodel_output_id => $genome_hash->{$genomes->[$i]}.".mdl",
					media_id => "Carbon-D-Glucose",
					template_id => $tid,
					template_workspace => "chenry:narrative_1493181437626",
					genome_workspace => "jplfaria:narrative_1492808527866",
					media_workspace => "KBaseMedia",
					gapfill_model => 1
				});
			}
		};
		if ($@) {
			my $error = $@;
			print $i.":".$genomes->[$i].":".$error."\n";
		} else {
			print $i.":".$genomes->[$i].":success\n";
		}
		$count++;
	}
}