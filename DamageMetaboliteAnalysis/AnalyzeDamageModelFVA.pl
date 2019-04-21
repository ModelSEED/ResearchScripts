use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "jplfaria:narrative_1510947016191";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

Bio::KBase::ObjectAPI::functions::func_run_flux_balance_analysis({
	workspace => $workspace,
	fbamodel_id => "CombinedModel",
	fba_output_id => "CombinedModel.fva",
	fbamodel_workspace => $workspace,
	media_id => "Complete",
	media_workspace => "KBaseMedia",
	minimize_flux => 1,
	fva => 1
});

my $fba = $impl->util_get_object(Bio::KBase::utilities::buildref("CombinedModel.fva",$workspace));

my $rxns = $fba->FBAReactionVariables();

print "ID\tMin\tMax\n";
for (my $i=0; $i < @{$rxns}; $i++) {
	if ($rxns->[$i]->max() > 0.000001 || $rxns->[$i]->min() < -0.000001) {
		my $mdlrxn = $rxns->[$i]->modelreaction();
		if ($mdlrxn->id() =~ /damage/) {
			print $mdlrxn->id()."\t".$rxns->[$i]->min()."\t".$rxns->[$i]->max()."\n";
		}
	}
}
for (my $i=0; $i < @{$rxns}; $i++) {
	if ($rxns->[$i]->max() > 0.000001 || $rxns->[$i]->min() < -0.000001) {
		my $mdlrxn = $rxns->[$i]->modelreaction();
		if ($mdlrxn->id() =~ /repair/) {
			print $mdlrxn->id()."\t".$rxns->[$i]->min()."\t".$rxns->[$i]->max()."\n";
		}
	}
}