use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/MetabolomicsData.tsv");
my $formula_hash = {};
my $peak_hash = {};
for (my $i=0; $i < @{$lines}; $i++) {
	my $array = [split(/[\t]/,$lines->[$i])];
	$peak_hash->{$array->[0]} = [$array->[1],$array->[2]];
	if (length($array->[2]) > 0) {
		$formula_hash->{$array->[2]}->{$array->[0]} = $array->[1];
	}
}

my $eligible_count = 0;
my $full_count = 0;
my $smiles_list = [];
my $smiles_hash = {};
my $compound_output = ["inchikey\tsmiles\tid\tsource"];
my $cpd_hash = Bio::KBase::utilities::compound_hash();
foreach my $cpd (keys(%{$cpd_hash})) {
	$full_count++;
	if (defined($cpd_hash->{$cpd}->{smiles})) {
		$eligible_count++;
		push(@{$smiles_list},$cpd_hash->{$cpd}->{smiles});
		$smiles_hash->{$cpd_hash->{$cpd}->{smiles}}->{$cpd} = 1;
		push(@{$compound_output},$cpd_hash->{$cpd}->{inchikey}."\t".$cpd_hash->{$cpd}->{smiles}."\t".$cpd_hash->{$cpd}->{id}."_c0\tDB");
	}
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/DB.compounds",$compound_output);
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/temp/smiles.smi",$smiles_list);
system('obabel /Users/chenry/temp/smiles.smi -oinchi -T /nochg/formula > /Users/chenry/temp/inchi');
$lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/temp/inchi");
my $hit_count = 0;
my $peak_count = 0;
my $peak_hits = {};
my $output = ["id\tmass\tformula\tmsids"];
for (my $i=0; $i < @{$lines}; $i++) {
	if ($lines->[$i] =~ m/InChI.1S\/([^\/]+)\//) {
		my $formula = $1;
		if (defined($formula_hash->{$formula})) {
			$hit_count++;
			foreach my $peakid (keys(%{$formula_hash->{$formula}})) {
				foreach my $cpdid (keys(%{$smiles_hash->{$smiles_list->[$i]}})) {
					$peak_hits->{$peakid}->{$cpdid."_c0"} = 1;
				}
			}
		}
	}
}
foreach my $peakid (keys(%{$peak_hits})) {
	$peak_count++;
	push(@{$output},$peakid."\t".$peak_hash->{$peakid}->[0]."\t".$peak_hash->{$peakid}->[1]."\t".join(";",keys(%{$peak_hits->{$peakid}})));	
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/DB.peaks",$output);

print "Peak count:".$peak_count."\n";
print "Compound count:".$hit_count."\n";
print "Eligible:".$eligible_count."\n";