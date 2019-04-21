use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("CombinedModel",$workspace));
my $cpds = $model->modelcompounds();

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/Peaklist.txt");
my $peakstring = "";
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	#First check for ModelSEED ID, then check for inchi-key
	if (defined($model->getObject("modelcompounds",$array->[1]))) {
		if (length($peakstring) > 0) {
			$peakstring .= ";";
		}
		$peakstring .= "peak.".$i.":".$array->[1];
	} else {
		#Then check if there is an inchi string match
		for (my $j=0; $j < @{$cpds}; $j++) {
			if ($cpds->[$j]->inchikey() eq $array->[0]) {
				if (length($peakstring) > 0) {
					$peakstring .= ";";
				}
				$peakstring .= "peak.".$i.":".$cpds->[$j]->id();
				last;
			}
		}
	}
}

print $peakstring."\n";

my $output = Bio::KBase::ObjectAPI::functions::func_gapfill_metabolic_model({	
	fbamodel_id => "iMB155.trans.noex.fix",
	workspace => $workspace,
	media_id => "Complete",
	media_workspace => "KBaseMedia",
	source_fbamodel_id => "CombinedModel",
	target_reaction => "bio1",
	fbamodel_output_id => "iMB155.trans.noex.fix.metgf",
	metabolomics_peak_data => $peakstring
});