#!/usr/bin/perl -w

use strict;
use warnings;
use Config::Simple;
use Bio::KBase::workspaceService::Client;
use Bio::KBase::fbaModelServices::Impl;

$|=1;
my $config = $ARGV[0];
my $overwrite = 0;
if (!defined($config)) {
	print STDERR "No config file provided!\n";
	exit(-1);
}
if (!-e $config) {
	print STDERR "Config file ".$config." not found!\n";
	exit(-1);
}
#Params: writesbml.wsurl, writesbml.fbaurl, writesbml.auth
my $c = Config::Simple->new();
$c->read($config);
my $wserv = Bio::KBase::workspaceService::Client->new($c->param("kbclientconfig.wsurl"));
my $fbaserv = Bio::KBase::fbaModelServices::Impl->new({workspace => $wserv});

my $directory = "/Users/chenry/workspace/PhenotypeData/";
my $genomes = [qw(
	62977.3 71421.1 83332.1 83333.1 85962.1 93062.4 158879.1 169963.1
	171101.1 208964.1 224308.1 243273.1 243277.1 266834.1 99287.1
	267608.1 269483.3 272626.1 272635.1 316407.3 345073.6 401614.5
)];

foreach my $genome (@{$genomes}) {
	my $dir = $directory.$genome."/";
	if (-e $dir."CultureConditions.txt") {
		print $genome." biolog\n";
		my $source;
		my $phenotypes = [];
		open(my $fh, "<".$dir."CultureConditions.txt") || die "could not open biolog file\n!"; 
		my $line = <$fh>;#Heading line
		while ($line = <$fh>) {
			chomp($line);
			my $items = [split(/\t/,$line)];
			if(!defined($source)) {
				$source = $items->[2];
			}
			my $media = $items->[0];
			$media =~ s/[\s,]/-/g;
			$media =~ s/[\'\(\)\`]//g;
			$media =~ s/--/-/g;
			if ($media !~ m/[\w-]/) {
				print "Bad media:".$genome.":".$media."\n";
			}
			push(@{$phenotypes},[[],$media,"KBaseMedia",[],$items->[1]]);
		}
		close($fh);
		if ($source =~ m/PMID(\d+)/) {
			$source = "PMID:".$1;
		}
		$fbaserv->import_phenotypes({
			phenotypeSet => $genome."-biolog",
			workspace => "KBasePhenotypeDatasets",
			genome => $genome,
			genome_workspace => "KBasePhenotypeDatasets",
			phenotypes => $phenotypes,
			type => "BinaryBiolog",
			name => $genome." biolog data",
			source => $source,
			auth => $c->param("kbclientconfig.auth"),
			ignore_errors => 1
		});
	}
	if (-e $dir."Essentiality.txt") {
		print $genome." essentiality\n";
		my $source;
		my $phenotypes = [];
		open(my $fh, "<".$dir."Essentiality.txt") || die "could not open biolog file\n!"; 
		my $line = <$fh>;#Heading line
		while ($line = <$fh>) {
			chomp($line);
			my $items = [split(/\t/,$line)];
			if(!defined($source)) {
				$source = $items->[3];
			}
			my $media = $items->[1];
			$media =~ s/[\s,]/-/g;
			$media =~ s/[\'\(\)\`]//g;
			if ($media !~ m/[\w-]/) {
				print "Bad media:".$media."\n";
			}
			if ($items->[2] =~ m/^essential/ || $items->[2] eq "potentially essential" || $items->[2] eq "potential essential") {
				push(@{$phenotypes},[[$items->[0]],$media,"KBaseMedia",[],0]);
			} elsif ($items->[2] =~ m/^nonessential/) {
				push(@{$phenotypes},[[$items->[0]],$media,"KBaseMedia",[],1]);
			} else {
				print "Unrecognized essentiality:".$genome.":".$items->[2]."\n";
			}
		}
		close($fh);
		if ($source =~ m/PMID(\d+)/) {
			$source = "PMID:".$1;
		}
		$fbaserv->import_phenotypes({
			phenotypeSet => $genome."-essentiality",
			workspace => "KBasePhenotypeDatasets",
			genome => $genome,
			genome_workspace => "KBasePhenotypeDatasets",
			phenotypes => $phenotypes,
			type => "BinaryGeneEssentiality",
			name => $genome." essentiality data",
			source => $source,
			auth => $c->param("kbclientconfig.auth"),
			ignore_errors => 1
		});
	}
}

1;
