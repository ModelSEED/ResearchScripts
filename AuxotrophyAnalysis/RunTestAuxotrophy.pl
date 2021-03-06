use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $output = Bio::KBase::kbaseenv::list_objects({
	workspaces => ["cshgusto:1484552981303"]
});

my $genomelist;
for (my $i=0; $i < @{$output}; $i++) {
	if ($output->[$i]->[2] =~ m/KBaseGenomes\.Genome/) {
		push(@{$genomelist},$output->[$i]->[1]);
	}
}

Bio::KBase::ObjectAPI::functions::func_predict_auxotrophy({
	workspace => "chenry:narrative_1508049563066",
	genome_ids => $genomelist,
	genome_workspace => "cshgusto:1484552981303"
});