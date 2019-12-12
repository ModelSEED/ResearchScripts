use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $filename = "/Users/chenry/Dropbox/workspace/PNNLSFA/Cjaponicus_fulldata_peaks.tsv";
my $lines = Bio::KBase::ObjectAPI::utilities::LOADFILE($filename);
my $output = ["id\tmass\tformula"];
for (my $i=1; $i < @{$lines}; $i++) {
	my $array = [split(/[\t]/,$lines->[$i])];
	my $id = "peak.".$i;
	my $mass = $array->[0];
	my $formula = $array->[2];
	if (defined($formula) && length($formula) > 0) {
		my $carbon = 0;
		if ($formula =~ m/13C/) {
			$carbon = 1;
			$formula =~ s/13C//;
		}
		if ($formula =~ m/(.*)C(\d+)([A-Z]*.*)/) {
			$carbon += $2;
			$formula = $1."C".$carbon.$3;
		} elsif ($formula =~ m/(.*)C([A-Z]*.*)/) {
			$carbon += 1;
			$formula = $1."C".$carbon.$3;
		} elsif ($carbon == 1) {
			$formula = "C".$formula;
		}
	}
	push(@{$output},$id."\t".$mass."\t".$formula);
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PNNLSFA/MetabolomicsData.tsv",$output);