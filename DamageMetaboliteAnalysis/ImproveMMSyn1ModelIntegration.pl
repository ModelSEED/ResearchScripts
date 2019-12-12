use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $model = $impl->util_get_object(Bio::KBase::utilities::buildref("MMSyn1.raw",50558));
my $input = {
	compartment => "c",
	compartment_index => 0,
	priority => 0,
	hashes => {
		ids => {},
		names => {},
		structures => {},
		base_structures => {},
		formulas => {}
	}
};
Bio::KBase::utilities::metabolite_hash($input);

#Replacing missing data
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/JCVI modeling/Published MM model/cpd_data.tsv");
print "ID\tCurrent match\tKEGG match\tInchimatch\n";
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	my $cpd = $model->getObject("modelcompounds","M_".$array->[0]);
	if (defined($cpd)) {
		my $cpdhash;
		my $sourcehash;
		$cpd->charge($array->[3]);
		$cpd->formula($array->[2]);
		$cpd->name($array->[1]);
		if (defined($cpd->dblinks()->{ModelSeed})) {
			for (my $j=0; $j < @{$cpd->dblinks()->{ModelSeed}}; $j++) {
				$cpdhash->{$cpd->dblinks()->{ModelSeed}->[$j]}->{c} = 1;
				$sourcehash->{c}->{$cpd->dblinks()->{ModelSeed}->[$j]} = 1;
			}
		}
		if ($array->[4] =~ m/C\d+/) {
			if (defined($input->{hashes}->{ids}->{$array->[4]})) {
				foreach my $newcpd (keys(%{$input->{hashes}->{ids}->{$array->[4]}})) {
					if ($newcpd =~ m/(cpd\d+)/) {
						$newcpd = $1;
						$cpdhash->{$newcpd}->{k} = 1;
						$sourcehash->{k}->{$newcpd} = 1;
					}
				}
			}
		}
		if (length($array->[5]) > 0) {
			$cpd->inchikey($array->[5]);
			if (defined($input->{hashes}->{structures}->{$array->[5]})) {
				foreach my $newcpd (keys(%{$input->{hashes}->{structures}->{$array->[5]}})) {
					if ($newcpd =~ m/(cpd\d+)/) {
						$newcpd = $1;
						$cpdhash->{$newcpd}->{s} = 1;
						$sourcehash->{s}->{$newcpd} = 1;
					}
				}
			}
		}
		print $array->[0]."\t";
		if (defined($cpd->dblinks()->{ModelSeed})) {
			my $output = "";
			for (my $j=0; $j < @{$cpd->dblinks()->{ModelSeed}}; $j++) {
				if (length($output) > 0) {
					$output .= "|";
				}
				$output .= $cpd->dblinks()->{ModelSeed}->[$j].":".join("",keys(%{$cpdhash->{$cpd->dblinks()->{ModelSeed}->[$j]}}));
			}
			print $output."\t";
		} else {
			print "none\t";
		}
		if ($array->[4] =~ m/C\d+/) {
			$cpd->dblinks()->{KEGG} = $array->[4];
			if (defined($input->{hashes}->{ids}->{$array->[4]})) {
				my $output = "";
				foreach my $newcpd (keys(%{$input->{hashes}->{ids}->{$array->[4]}})) {
					if ($newcpd =~ m/(cpd\d+)/) {
						$newcpd = $1;
						if (length($output) > 0) {
							$output .= "|";
						}
						$output .= $newcpd.":".join("",keys(%{$cpdhash->{$newcpd}}));	
					}
				}
				print $output."\t";
			} else {
				print "not in db\t";
			}
		} else {
			print "none\t";
		}
		if (length($array->[5]) > 0) {
			if (defined($input->{hashes}->{structures}->{$array->[5]})) {
				my $output = "";
				foreach my $newcpd (keys(%{$input->{hashes}->{structures}->{$array->[5]}})) {
					if ($newcpd =~ m/(cpd\d+)/) {
						$newcpd = $1;
						if (length($output) > 0) {
							$output .= "|";
						}
						$output .= $newcpd.":".join("",keys(%{$cpdhash->{$newcpd}}));	
					}
				}
				print $output."\n";
			} else {
				print "not in db\n";
			}
		} else {
			print "none\n";
		}
	}
}
exit;
#Fixing IDs
my $cpds = $model->modelcompounds();
my $rxns = $model->modelreactions();
my $translation = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->id();
	$id =~ s/^M_//g;
	$id =~ s/[^\w]/-/g;
	$id =~ s/\s/-/g;
	$id .= "0";
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
my $wsmeta = $impl->util_save_object($model,"50558/MMSyn1.fix");