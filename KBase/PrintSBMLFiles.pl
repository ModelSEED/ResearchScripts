use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;
my $directory = $ARGV[0];

Bio::KBase::kbaseenv::create_context_from_client_config();
my $ws = Bio::KBase::kbaseenv::ws_client();
my $impl = fba_tools::fba_toolsImpl->new();
my $store = $impl->util_store();
my $models = $ws->list_objects({
	workspaces => ["janakakbase:CM-VR"],
	type => "KBaseFBA.FBAModel",
});
for (my $i=0; $i < @{$models}; $i++) {
	my $model = $store->get_object("janakakbase:CM-VR/".$models->[$i]->[1]);
	my $sbml = $model->export({format => "sbml"});
	open ( my $fh, ">", $directory."/".$models->[$i]->[1]);
    print $fh $sbml."\n";
    close($fh);
}