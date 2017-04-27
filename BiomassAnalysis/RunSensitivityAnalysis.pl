use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $procs = $ARGV[0];
my $index = $ARGV[1];

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();

my $modelws = "20905";
my $fbaws = "20955";
my $media = "Carbon-D-Glucose";
my $mediaws = "KBaseMedia";

my $ws = Bio::KBase::kbaseenv::ws_client();
my $models = $ws->list_objects({
	workspace => [$modelws],
	type => "KBaseFBA.FBAModel"
});

for (my $i=0; $i < @{$models}; $i++) {
	if ($i % $procs  == $index) {
		Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
			workspace => $fbaws,
			fbamodel_id => $models->[$i]->[1],
			fba_output_id => $models->[$i]->[1].".sensfba",
			fbamodel_workspace => $modelws,
			media_id => $media,
			media_workspace => $mediaws,
			minimize_flux => 1,
			sensitivity_analysis => 1
		});
		Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
			workspace => $fbaws,
			fbamodel_id => $models->[$i]->[1],
			fba_output_id => $models->[$i]->[1].".mmfva",
			fbamodel_workspace => $modelws,
			media_id => $media,
			media_workspace => $mediaws,
			minimize_flux => 1,
			fva => 1
		});
		Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
			workspace => $fbaws,
			fbamodel_id => $models->[$i]->[1],
			fba_output_id => $models->[$i]->[1].".comfva",
			fbamodel_workspace => $modelws,
			minimize_flux => 1,
			fva => 1
		});
	}
}