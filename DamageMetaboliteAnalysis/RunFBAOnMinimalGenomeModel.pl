use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.trans.noex.fix",$workspace));

Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
	workspace => $workspace,
	fbamodel_id => "iMB155.trans.noex.fix",
	fba_output_id => "iMB155.trans.noex.fix.fba"
})