#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(getUser printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

if (!defined($ARGV[2])) {
	print "Insufficient arguments provided. Usage:\nperl load_contigs.pl <fasta file> <workspace id> <contig set name>\n";
}

my $fastafile = $ARGV[0];
my $workspace = $ARGV[1];
my $contig_name = $ARGV[2];

my $object = {
	id => $contig_name,
	name => $contig_name,
	source_id => $contig_name,
	source => getUser(),
	type => "Organism",
	contigs => [],
	md5 => undef
};
open FILE, "<",$fastafile  or do {
	die "$0: open ".$fastafile.": $!";
};
my $contig;
my $sequence = "";
my $description;
my $length = 0;
my $lengths = [];
my $contigids = [];
my $contigseq = {};
my $gc = 0;
my $numcontigs = 0;
while (my $line = <FILE>) {
	chomp($line);
	if ($line =~ m/\>(.+)/) {
		if (defined($contig)) {
			$numcontigs++;
			push(@{$lengths},length($sequence));
			push(@{$contigids},$contig);
			$length += length($sequence);
			$contigseq->{$contig} = $sequence;
			push(@{$object->{contigs}},{
				id => $contig,
				"length" => length($sequence),
				md5 => Digest::MD5::md5_hex($sequence),
				sequence => $sequence,
				genetic_code => 11,
				name => $contig,
				complete => 1,
				description => $description
			});
			$sequence =~ s/[atAT]//g;
			$gc += length($sequence);
		}
		$sequence = "";
		$contig = $1;
		$description = $1;
		if($contig =~ /(.*)?\s(.+)/ ) {
			$contig = $1;
			$description = $2;
		}
	} else {
		$sequence .= $line;
	}
}
if (defined($contig)) {
	$numcontigs++;
	push(@{$lengths},length($sequence));
	push(@{$contigids},$contig);
	$length += length($sequence);
	$contigseq->{$contig} = $sequence;
	push(@{$object->{contigs}},{
		id => $contig,
		"length" => length($sequence),
		md5 => Digest::MD5::md5_hex($sequence),
		sequence => $sequence,
		genetic_code => 11,
		name => $contig,
		complete => 1,
		description => $description
	});
	$sequence =~ s/[atAT]//g;
	$gc += length($sequence);
}
$gc = $gc/$length;
close(FILE);
$object->{contigs} = [sort { $a->{sequence} cmp $b->{sequence} } @{$object->{contigs}}];
my $str = "";
for (my $i=0; $i < @{$object->{contigs}}; $i++) {
	if (length($str) > 0) {
		$str .= ";";
	}
	$str .= $object->{contigs}->[$i]->{sequence};
}
$object->{md5} = Digest::MD5::md5_hex($str);
save_workspace_object($workspace."/".$object->{id},$object,"KBaseGenomes.ContigSet");