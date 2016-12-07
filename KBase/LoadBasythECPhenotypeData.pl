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
my $phenotypes = [];
open(my $fh, "</Users/chenry/workspace/BaSynthEC/Phenotypes.txt"); 
my $line = <$fh>;#Heading line
while ($line = <$fh>) {
	chomp($line);
	my $items = [split(/\t/,$line)];
	push(@{$phenotypes},[[split(/;/,$items->[3])],$items->[1],"BaSynthECPhenotypes",[],$items->[2],$items->[0]]);
}
close($fh);
$fbaserv->import_phenotypes({
	workspace => "BaSynthECPhenotypes",
	genome => "BSubtilis168",
	genome_workspace => "BaSynthECGenome",
	phenotypes => $phenotypes,
	type => "BinaryMultiKO",
	name => "BaSynthEC interval phenotypes",
	source => "PMID:23109554",
	auth => $c->param("kbclientconfig.auth"),
	ignore_errors => 1
});

1;
