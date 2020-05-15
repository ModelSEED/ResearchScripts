use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $workspace = "chenry:narrative_1520492239994";

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($impl);

my $atomhash = Bio::KBase::constants::atomic_masses();
my $cpdhash = Bio::KBase::utilities::compound_hash();

my $template = $impl->util_get_object("NewKBaseModelTemplates/GramNegModelTemplateV2");
my $model = $impl->util_get_object("46377/4/1");
my $comps = $model->modelcompounds();
my $modelhash = {};
for (my $i=0; $ i < @{$comps}; $i++) {
	$modelhash->{$comps->[$i]->msid()} = $comps->[$i]->neutral_formula();
}
my $templatehash = {};
$comps = $template->compcompounds();
for (my $i=0; $ i < @{$comps}; $i++) {
	if ($comps->[$i]->id() =~ m/(cpd\d+)/ && defined($comps->[$i]->formula())) {
		$templatehash->{$1} = $comps->[$i]->neutral_formula();
	}
}

print "Peak name\tSEED\tFormula\tCalc mass\tNeutral mass\tM/Z\t\t\tCalculate mass\tSEED formula\tSEED mass\tSEED calc mass\tSEED charge\tTemp NF\tModel NF\tSEED NF\n";
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE("/Users/chenry/Dropbox/workspace/Metabolite repair/WOMPeaks.txt");
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/\t/,$lines->[$i])];
	my $seedids = [split(/;/,$array->[3])];
	for (my $j=0; $j < @{$seedids}; $j++) {
		my $cpddata = $cpdhash->{$seedids->[$j]};
		my $calc_seed_mass = 0;
		my $calc_in_formula = 0;
		my $atoms = ["C","N","O","P","S","H"];
		for (my $k=0; $k < @{$atoms}; $k++) {
			my $a = $atoms->[$k];
			my $count = 0;
			if ($array->[1] =~ m/$a(\d+)/) {
				$count = $1;
			} elsif ($array->[1] =~ m/$a/) {
				$count = 1;
			}
			$calc_in_formula += $atomhash->{$a}*$count;
			$count = 0;
			if ($cpddata->{formula} =~ m/$a(\d+)/) {
				$count = $1;
			} elsif ($cpddata->{formula} =~ m/$a/) {
				$count = 1;
			}
			$calc_seed_mass += $atomhash->{$a}*$count;
		}
		print $array->[0]."\t".$seedids->[$j]."\t".$array->[1]."\t".$array->[2]."\t".$array->[4]."\t".$array->[5]."\t".$array->[6]."\t".$array->[7]."\t".$calc_in_formula."\t".$cpddata->{formula}."\t".$cpddata->{mass}."\t".$calc_seed_mass."\t".$cpddata->{charge}."\t";
		if (defined($templatehash->{$seedids->[$j]})) {
			print $templatehash->{$seedids->[$j]};
		} else {
			print "";
		}
		print "\t";
		if (defined($modelhash->{$seedids->[$j]})) {
			print $modelhash->{$seedids->[$j]};
		} else {
			print "";
		}
		print "\t".$cpddata->{neutral_formula}."\n";
	} 
}