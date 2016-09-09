use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;

my $procs = $ARGV[0];
my $index = $ARGV[1];

my $impl = fba_tools::fba_toolsImpl->new();
my $ws = $impl->util_kbase_store()->workspace();
my $genomes = $ws->list_objects({
	workspaces => ["coremodels"],
	type => "KBaseGenomes.Genome",
});

my $count = 0;
for (my $i=0; $i < @{$genomes}; $i++) {
	if ($i % $procs  == $index) {
		if ($count % 100 == 0) {
			$impl = fba_tools::fba_toolsImpl->new();		}
		eval {
			$impl->build_metabolic_model({
				workspace => "PubSEEDGramNegModels",
				genome_id => $genomes->[$i]->[1],
				fbamodel_output_id => $genomes->[$i]->[1].".model",
				media_id => "Carbon-D-Glucose",
				template_id => "gramneg",
				genome_workspace => "coremodels",
				media_workspace => "KBaseMedia",
				gapfill_model => 1
			});
		};
		if ($@) {
			my $error = $@;
			print $i.":".$genomes->[$i]->[1].":".$error;
		} else {
			print $i.":".$genomes->[$i]->[1].":success";
		}
		$count++;
	}
}