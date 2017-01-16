use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;
my $directory = $ARGV[0];

Bio::KBase::kbaseenv::create_context_from_client_config();
my $ws = Bio::KBase::kbaseenv::ws_client();
my $impl = fba_tools::fba_toolsImpl->new();
my $store = $impl->util_store();
my $genomes = $ws->list_objects({
	workspaces => ["PamGenomes"],
	type => "KBaseGenomes.Genome",
});
for (my $i=0; $i < @{$genomes}; $i++) {
	my $genome = $store->get_object($genomes->[$i]->[7]."/".$genomes->[$i]->[1]);
	my $completeness = $genome->compute_genome_completeness();
	print $genomes->[$i]->[1]."\t".$completeness."\n";
}