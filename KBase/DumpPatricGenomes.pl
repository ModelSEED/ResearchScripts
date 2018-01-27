#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use P3DataAPI;
use gjoseqlib;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $filename = $ARGV[0];

my $output = {};

my $d = P3DataAPI->new();
$d->{chunk_size} = 5000;
my @res = $d->query("genome", ["ne", "genome_id", "test"],
		    ["select", "genome_id", "genome_name", "taxon_id", "taxon_lineage_names"],
		   );

my $count = @res;
print "Total genomes = ".$count."\n";

for (my $i=0; $i < $count; $i++) {
	if (defined($res[$i])) {
		my $tax = $res[$i]->{taxon_lineage_names};
		shift @$tax if $tax->[0] eq 'root' || $tax->[0] eq 'cellular organisms';
		my $dom = $tax->[0];
		my $domain = $dom =~ /^Arch/i      ? 'Archaea'
		    : $dom =~ /^Bact/i      ? 'Bacteria'
		    : $dom =~ /^Eu[ck]a/i   ? 'Eukaryota'
		    : $dom =~ /^Vir/i       ? 'Virus'
		    : $dom =~ /^Environ/i   ? 'Environmental Sample'
		    : $dom =~ /^Un\S* Env/i ? 'Environmental Sample'
		    :                         'Unknown';
		my $taxonomy = join("; ", @$tax, $res[$i]->{genome_name});
		$output->{$res[$i]->{genome_id}} = {
			ti => $res[$i]->{taxon_id},
			n => $res[$i]->{genome_name},
			d => $dom,
			t => $taxonomy,
			ssu => []
		};
	}
}

print "Now retrieving 16s!\n";

my(@fields) = qw(feature_type
		 patric_id
		 na_sequence
		 product
		 );
		 
@res = $d->query("genome_feature",
	["select", @fields],
	["eq", "feature_type", "rRNA"],
	["eq", "annotation", "PATRIC"],
	["eq", "product", "SSU rRNA"],
);

my $count = @res;
print "Total 16s = ".$count."\n";

for (my $i=0; $i < $count; $i++) {
	if ($res[$i]->{patric_id} =~ m/fig\|(\d+\.\d+)\./) {
		my $genomeid = $1;
		if (defined($output->{$genomeid})) {
			push(@{$output->{$genomeid}->{ssu}},$res[$i]->{na_sequence});
		}
	}
}

Bio::KBase::ObjectAPI::utilities::PRINTFILE($filename,[Bio::KBase::ObjectAPI::utilities::TOJSON($output,1)]);