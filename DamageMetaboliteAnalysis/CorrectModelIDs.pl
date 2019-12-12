use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn1.raw",$workspace));

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/Published MM model/cpd_data.tsv");
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	my $cpd = $model->getObject("modelcompounds","M_".$array->[0]);
	$cpd->id()
	
	my $kegg = $array->[1];
	my $charge
	
	
	
	$cpd->inchikey($array->[1]);
}





my $cpds = $model->modelcompounds();
my $rxns = $model->modelreactions();
my $translation = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->id();
	$id =~ s/^M_//g;
	$id =~ s/[^\w]/-/g;
	$id =~ s/\s/-/g;
	$translation->{$cpds->[$i]->id()} = $id;
	$cpds->[$i]->id($id)
}
my $bios = $model->biomasses();
for (my $i=0; $i < @{$bios}; $i++) {
	my $biocpds = $bios->[$i]->biomasscompounds();
	for (my $j=0; $j < @{$biocpds}; $j++) {
		my $biocpd = $biocpds->[$j];
		if ($biocpd->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$biocpd->modelcompound_ref($1.$translation->{$2});
			}
		}
	}
}
for (my $i=0; $i < @{$rxns}; $i++) {
	$rxns->[$i]->maxforflux(0);
	$rxns->[$i]->maxrevflux(0);
	if ($rxns->[$i]->direction() eq ">" || $rxns->[$i]->direction() eq "=") {
		$rxns->[$i]->maxforflux(1000);
	}
	if ($rxns->[$i]->direction() eq "<" || $rxns->[$i]->direction() eq "=") {
		$rxns->[$i]->maxrevflux(1000);
	}
	my $rgts = $rxns->[$i]->modelReactionReagents();
	for (my $j=0; $j < @{$rgts}; $j++) {
		if ($rgts->[$j]->modelcompound_ref() =~ m/(.+\/)([^\/]+$)/) {
			if (defined($translation->{$2})) {
				$rgts->[$j]->modelcompound_ref($1.$translation->{$2});
			}
		}
	}	
}

my $wsmeta = $impl->util_save_object($model,$workspace."/MMSyn1.fix");