use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn1.int",50558));
my $cpdhash = Bio::KBase::utilities::compound_hash();
my $translation = {};

#Replacing missing data
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/Published MM model/cpd_data.tsv");
my $inchihash = {};
my $count = 0;
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	my $cpd = $model->getObject("modelcompounds","M_".$array->[0]);
	if (defined($cpd)) {
		$count++;
		$cpd->charge($array->[3]);
		$cpd->formula($array->[2]);
		$cpd->name($array->[1]);
		if (length($array->[6]) > 0) {
			$translation->{$cpd->id()} = $array->[6]."_".$cpd->modelcompartment()->id();
			$cpd->id($array->[6]."_".$cpd->modelcompartment()->id());
			$cpd->dblinks()->{ModelSeed} = [$array->[6]];
			if (defined($cpdhash->{$array->[6]}->{smiles})) {
				$cpd->smiles($cpdhash->{$array->[6]}->{smiles});
			}
			if (defined($cpdhash->{$array->[6]}->{inchikey})) {
				#$inchihash->{$cpdhash->{$array->[6]}->{inchikey}} = $cpd;
				$cpd->inchikey($cpdhash->{$array->[6]}->{inchikey});
			}
		}
		if ($array->[4] =~ m/C\d+/) {
			$cpd->dblinks()->{KEGG} = [$array->[4]];
		}
		if (length($array->[5]) > 25 && length($cpd->inchikey()) == 0) {
			$inchihash->{$array->[5]} = $cpd;
			$cpd->inchikey($array->[5]);
		}
	} else {
		print "Fail:M_".$array->[0]."\n";
	}
}
print "Total:".$count."\n";

#Querying pubchem to get smiles from inchikey
print "Inchi:\n".join("\n",keys(%{$inchihash}))."\n";
my $output = Bio::KBase::ObjectAPI::utilities::query_pubchem({inchikey => [keys(%{$inchihash})]});
print Bio::KBase::ObjectAPI::utilities::TOJSON($output);
foreach my $inchikey (keys(%{$inchihash})) {
	if (defined($output->{$inchikey})) {
		$inchihash->{$inchikey}->smiles($output->{$inchikey}->{CanonicalSMILES});
	}
}

#Fixing IDs
my $cpds = $model->modelcompounds();
my $rxns = $model->modelreactions();
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->id();
	if ($id !~ m/cpd\d+/) {
		$id =~ s/^M_//g;
		$id =~ s/[^\w]/-/g;
		$id =~ s/\s/-/g;
		$id =~ s/_([a-z])$/=$1/;
		$id =~ s/_//g;
		$id =~ s/=([a-z])$/_$1/;
		$id .= "0";
		$translation->{$cpds->[$i]->id()} = $id;
		$cpds->[$i]->id($id);
	}
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
	if ($rxns->[$i]->reaction_ref() =~ m/(rxn\d+)/) {
		my $rxnid = $1;
		$rxns->[$i]->reaction_ref("~/template/reactions/id/".$rxnid."_c");
		if ($rxnid ne "rxn00000") {
			$rxns->[$i]->id($rxnid."_c0");
		}
	}
	if ($rxns->[$i]->id() !~ m/rxn\d+/) {	
		my $id = $rxns->[$i]->id();
		$id =~ s/^R_//g;
		$id =~ s/[^\w]/-/g;
		$id =~ s/\s/-/g;
		$id =~ s/_/-/g;
		$id .= "_c0";
		$rxns->[$i]->id($id);
	}
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
my $wsmeta = $impl->util_save_object($model,"50558/MMSyn1");