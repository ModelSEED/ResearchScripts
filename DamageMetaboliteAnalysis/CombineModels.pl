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
my $cpdhash = Bio::KBase::utilities::compound_hash();

my $cpd_id_hash;
my $inchi_key_hash;
my $cpds = $model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if ($cpds->[$i]->id() =~ m/(cpd\d+)/) {
		$cpd_id_hash->{$1} = "base";
	}
	if (length($cpds->[$i]->inchikey()) > 0) { 
		$inchi_key_hash->{$cpds->[$i]->inchikey()} = "base";
	} elsif ($cpds->[$i]->id() =~ m/(cpd\d+)/) {
		if (defined($cpdhash->{$1}->{inchikey})) {
			$inchi_key_hash->{$cpdhash->{$1}->{inchikey}} = "base";
		}
	}
}
$cpds = $dam_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if ($cpds->[$i]->id() =~ m/(cpd\d+)/ && !defined($cpd_id_hash->{$1})) {
		$cpd_id_hash->{$1} = "damage";
	}
	if (length($cpds->[$i]->inchikey()) > 0 && !defined($inchi_key_hash->{$cpds->[$i]->inchikey()})) { 
		$inchi_key_hash->{$cpds->[$i]->inchikey()} = "damage";
	} elsif ($cpds->[$i]->id() =~ m/(cpd\d+)/ && !defined($cpdhash->{$1})) {
		if (defined($cpdhash->{$1}->{inchikey}) && !defined($inchi_key_hash->{$cpdhash->{$1}->{inchikey}})) {
			$inchi_key_hash->{$cpdhash->{$1}->{inchikey}} = "damage";
		}
	}
}
$cpds = $rep_model->modelcompounds();
for (my $i=0; $i < @{$cpds}; $i++) {
	if ($cpds->[$i]->id() =~ m/(cpd\d+)/ && !defined($cpd_id_hash->{$1})) {
		$cpd_id_hash->{$1} = "repair";
	}
	if (length($cpds->[$i]->inchikey()) > 0 && !defined($inchi_key_hash->{$cpds->[$i]->inchikey()})) { 
		$inchi_key_hash->{$cpds->[$i]->inchikey()} = "repair";
	} elsif ($cpds->[$i]->id() =~ m/(cpd\d+)/ && !defined($cpdhash->{$1})) {
		if (defined($cpdhash->{$1}->{inchikey}) && !defined($inchi_key_hash->{$cpdhash->{$1}->{inchikey}})) {
			$inchi_key_hash->{$cpdhash->{$1}->{inchikey}} = "repair";
		}
	}
}

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI-metabolomics.txt");
print "Metabolite name	INCHIKEY	Chemical Name	Compound ID	InChI	Model\n";
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	print $lines->[$i]."\t";
	if (defined($cpd_id_hash->{$array->[3]})) {
		print $cpd_id_hash->{$array->[3]};
	} elsif (defined($inchi_key_hash->{$array->[1]})) {
		print $inchi_key_hash->{$array->[1]};
	}
	print "\n";
}