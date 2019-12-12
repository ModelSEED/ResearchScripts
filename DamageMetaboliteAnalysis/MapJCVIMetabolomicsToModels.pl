use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.trans.noex.fix",$workspace));
my $dam_model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.damage",$workspace));
my $rep_model = $impl->util_get_object(Bio::KBase::utilities::buildref("iMB155.repair",$workspace));
my $combo_model = $impl->util_get_object(Bio::KBase::utilities::buildref("CombinedModel",$workspace));
my $cpdhash = Bio::KBase::utilities::compound_hash();

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI-confident");
my $pubchem_args = {inchikey => []};
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	push	(@{$pubchem_args->{inchikey}},$array->[4]);
}
my $propertyhash = Bio::KBase::ObjectAPI::utilities::query_pubchem($pubchem_args);
my $full_hash = {};
my $db_base_inchi_hash = {};
foreach my $cpd (keys(%{$cpdhash})) {
	if (defined($cpdhash->{$cpd}->{inchikey})) {
		my $array = [split(/[_-]/,$cpdhash->{$cpd}->{inchikey})];
		$full_hash->{$array->[0]}->{$cpd} = 1;
		$db_base_inchi_hash->{$array->[0]}->{$cpd} = 1;
	}	
}

my $mdl_base_inchi_hash;
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (length($cpds->[$i]->inchikey()) > 0) { 
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		my $id = $cpds->[$i]->id();
		$id =~ s/_c0$//;
		$id =~ s/_e0$//;
		$full_hash->{$array->[0]}->{$id} = 1;
		$mdl_base_inchi_hash->{$array->[0]}->{$id} = 1;
	}
}

my $dam_base_inchi_hash;
$cpds = $dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (length($cpds->[$i]->inchikey()) > 0) { 
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		my $id = $cpds->[$i]->id();
		$id =~ s/_c0$//;
		$id =~ s/_e0$//;
		$full_hash->{$array->[0]}->{$id} = 1;
		$dam_base_inchi_hash->{$array->[0]}->{$id} = 1;
	}
}

my $rep_base_inchi_hash;
$cpds = $rep_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (length($cpds->[$i]->inchikey()) > 0) { 
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		my $id = $cpds->[$i]->id();
		$id =~ s/_c0$//;
		$id =~ s/_e0$//;
		$full_hash->{$array->[0]}->{$id} = 1;
		$rep_base_inchi_hash->{$array->[0]}->{$id} = 1;
	}
}

my $combo_base_inchi_hash;
$cpds = $combo_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if (length($cpds->[$i]->inchikey()) > 0) { 
		my $array = [split(/[_-]/,$cpds->[$i]->inchikey())];
		my $id = $cpds->[$i]->id();
		$id =~ s/_c0$//;
		$id =~ s/_e0$//;
		$full_hash->{$array->[0]}->{$id} = 1;
		$combo_base_inchi_hash->{$array->[0]}->{$id} = 1;
	}
}
print "ID\tName\tFormula\tInchikey\tMSI\tAdduct\tPlatform\tRT\tMZ\tSEED\tModel\tDamage\tRepair\tSEEDDB\n";
my $peaks = {};
my $count = 0;
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	my $inchi_array = [split(/[_-]/,$array->[4])];
	my $formula = "";
	if (defined($propertyhash->{$array->[4]}->{MolecularFormula})) {
		$formula = $propertyhash->{$array->[4]}->{MolecularFormula};
	}
	print $array->[0]."\t".$array->[1]."\t".$formula."\t".$array->[4]."\t".$array->[2]."\t".$array->[3]."\t".$array->[5]."\t".$array->[7]."\t".$array->[6]."\t";
	if (defined($full_hash->{$inchi_array->[0]})) {
		print join("|",keys(%{$full_hash->{$inchi_array->[0]}}))."\t"
	} else {
		print "\t";
	}
	if (defined($mdl_base_inchi_hash->{$inchi_array->[0]})) {
		print join("|",keys(%{$mdl_base_inchi_hash->{$inchi_array->[0]}}))."\t"
	} else {
		print "\t";
	}
	if (defined($dam_base_inchi_hash->{$inchi_array->[0]})) {
		print join("|",keys(%{$dam_base_inchi_hash->{$inchi_array->[0]}}))."\t"
	} else {
		print "\t";
	}
	if (defined($rep_base_inchi_hash->{$inchi_array->[0]})) {
		print join("|",keys(%{$rep_base_inchi_hash->{$inchi_array->[0]}}))."\t"
	} else {
		print "\t";
	}
	if (defined($db_base_inchi_hash->{$inchi_array->[0]})) {
		print join("|",keys(%{$db_base_inchi_hash->{$inchi_array->[0]}}))."\n"
	} else {
		print "\n";
	}
}