#!/usr/bin/perl -w

use strict;
use Config::Simple;
use Bio::KBase::fbaModelServices::Impl;
use File::Path;
$|=1;

my $config = $ARGV[0];
if (!defined($config)) {
	print STDERR "No config file provided!\n";
	exit(-1);
}
if (!-e $config) {
	print STDERR "Config file ".$config." not found!\n";
	exit(-1);
}

my $c = Config::Simple->new();
$c->read($config);

$Bio::KBase::fbaModelServices::Server::CallContext = {token => $c->param("kbclientconfig.auth")};
my $fba = Bio::KBase::fbaModelServices::Impl->new({"workspace-url" => "http://kbase.us/services/ws"});
$fba->_setContext($Bio::KBase::fbaModelServices::Server::CallContext,{});
my $ws = $fba->_workspaceServices();
my $genomes = $ws->list_objects({
	workspaces => ["KBasePublicGenomesV5"],
	type => "KBaseGenomes.Genome",
});

my $genomehash;
for (my $i=0; $i < @{$models}; $i++) {
	if ($genomes->[$i]->[1] =~ m/^(\d+\.\d+)\./) {
		$genomehash->{$1} = 1;
	}
}
for (my $i=0; $i < @{$genomes}; $i++) {
	if (!defined($genomehash->{$genomes->[$i]->[1]})) {
		$fba->queue_job({method => "genome_to_fbamodel",parameters => {
			templatemodel => "FullBiomassTemplate",
	    	templatemodel_workspace => "chenrydemo",
			genome => $genomes->[$i]->[1],
			workspace => "chenry:BiomassAnalysisMMModels",
	    	genome_workspace => "pubSEEDGenomes",
	    	model => $genomes->[$i]->[1].".fbamdl",
		}});
	}
};
