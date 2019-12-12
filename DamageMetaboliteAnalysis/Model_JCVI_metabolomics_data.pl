use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1566275241019";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $choice = $ARGV[0];

if ($choice == 1) {
	my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
		fbamodel_id => "iMB155.trans.noex.fix",
		fbamodel_workspace => "chenry:narrative_1520492239994",
		workspace => $workspace,
		media_id => "Complete",
		media_workspace => "KBaseMedia",
		source_fbamodel_id => "CombinedModel",
		source_fbamodel_workspace => "chenry:narrative_1520492239994",
		target_reaction => "bio1",
		fbamodel_output_id => "iMB155.combofit",
		metabolite_ref => "/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/NewDataset2.txt",
		metabolite_workspace => $workspace,
		metabolite_condition => "jcvi_syn3",
	});
} elsif ($choice == 2) {
	my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
		fbamodel_id => "iMB155.trans.noex.fix",
		fbamodel_workspace => "chenry:narrative_1520492239994",
		workspace => $workspace,
		media_id => "Complete",
		media_workspace => "KBaseMedia",
		target_reaction => "bio1",
		fbamodel_output_id => "iMB155.dbfix",
		metabolite_ref => "/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/NewDataset2.txt",
		metabolite_workspace => $workspace,
		metabolite_condition => "jcvi_syn3",
	});
} elsif ($choice == 3) {
	my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
		fbamodel_id => "iMB155.trans.noex.fix",
		fbamodel_workspace => "chenry:narrative_1520492239994",
		source_fbamodel_id => "CombinedModel",
		source_fbamodel_workspace => "chenry:narrative_1520492239994",
		workspace => $workspace,
		media_id => "Complete",
		media_workspace => "KBaseMedia",
		target_reaction => "bio1",
		fbamodel_output_id => "iMB155.chemofix",
		metabolite_ref => "/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/NewDataset2.txt",
		metabolite_workspace => $workspace,
		metabolite_condition => "jcvi_syn3",
		add_external_reactions => 0
	});
}