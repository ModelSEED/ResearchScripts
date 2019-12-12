#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
use Bio::KBase::utilities;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $workspace = $ARGV[0];
my $fba_id = $ARGV[1];
my $genome = $ARGV[2];

my $fba = $handler->util_get_object(Bio::KBase::utilities::buildref($fba_id,$workspace));
my $rxns = $fba->FBAReactionVariables();
my $db_hash = Bio::KBase::ObjectAPI::utilities::reaction_hash();

my $reaction_list = ["Reaction ID\tEC number\tReaction Eq\tReaction name\tReaction GPR\tFlux"];
my $gene_hash = {};
my $rxn_hash = {};
for (my $i=0; $i < @{$rxns}; $i++) {
	my $mdlrxn = $rxns->[$i]->modelreaction();
	$rxn_hash->{$mdlrxn->id()} = {
		equation => $mdlrxn->definition(),
		flux => $rxns->[$i]->value(),
		name => $mdlrxn->name(),
		gpr => $mdlrxn->gprString(),
		ec => {}
	};
	if ($mdlrxn->id() =~ m/(rxn\d+)/) {
		my $id = $1;
		if (defined($db_hash->{$mdlrxn->id()}->{ec_numbers})) {
			$rxn_hash->{$mdlrxn->id()}->{ec} = $db_hash->{$mdlrxn->id()}->{ec_numbers}->[0];
			$rxn_hash->{$mdlrxn->id()}->{ec} =~ s/^EC-//;
		} elsif (defined($db_hash->{$mdlrxn->id()}->{roles})) {
			for (my $j=0; $j < @{$db_hash->{$mdlrxn->id()}->{roles}}; $j++) {
				if ($db_hash->{$mdlrxn->id()}->{roles}->[$j] =~ m/([\d-]+\.[\d-]+\.[\d-]+\.[\d-]+)/) {
					$rxn_hash->{$mdlrxn->id()}->{ec}->{$1} = 1;
					last;
				}
			}
		}
	}
	push(@{$reaction_list},$mdlrxn->id()."\t".$rxn_hash->{$mdlrxn->id()}->{ec}."\t".$mdlrxn->definition()."\t".$mdlrxn->name()."\t".$mdlrxn->gprString()."\t".$rxns->[$i]->value());
	my $prots = $mdlrxn->modelReactionProteins();
	for (my $j=0; $j < @{$prots}; $j++) {
		my $subunits = $prots->[$j]->modelReactionProteinSubunits();
		for (my $k=0; $k < @{$subunits}; $k++) {
			my $ftrs = $subunits->[$k]->features();
			for (my $m=0; $m < @{$ftrs}; $m++) {
				if (!defined($gene_hash->{$ftrs->[$m]->id()})) {
					$gene_hash->{$ftrs->[$m]->id()} = {
						id => $ftrs->[$m]->id(),
						function => $ftrs->[$m]->function(),
						reactions => {}
					};
					if (length($gene_hash->{$ftrs->[$m]->id()}->{function}) == 0 && defined($ftrs->[$m]->functions()->[0])) {
						$gene_hash->{$ftrs->[$m]->id()}->{function} = join(" @ ",@{$ftrs->[$m]->functions()});
					} 
				}
				$gene_hash->{$ftrs->[$m]->id()}->{reactions}->{$mdlrxn->id()} = 1;
			}
		}
	}
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PlantSEED/CCRStudy/".$genome."_reaction_flux",$reaction_list);

my $gene_list = ["Gene ID\tFunction\tReactions\tHigh flux reaction\tEC number\tReaction Eq\tReaction name\tReaction GPR\tHighest flux\tTotal flux"];
foreach my $gene (keys(%{$gene_hash})) {
	my $line = $gene_hash->{$gene}->{id}."\t".$gene_hash->{$gene}->{function}."\t".join(";",keys(%{$gene_hash->{$gene}->{reactions}}))."\t";
	my $totalflux = 0;
	my $bestrxn;
	foreach my $rxnid (keys(%{$gene_hash->{$gene}->{reactions}})) {
		$totalflux += abs($rxn_hash->{$rxnid}->{flux});
		if (!defined($bestrxn) || abs($rxn_hash->{$rxnid}->{flux}) > abs($rxn_hash->{$bestrxn}->{flux})) {
			$bestrxn = $rxnid;
		}
	}
	$line .= $bestrxn."\t".$rxn_hash->{$bestrxn}->{ec}."\t".$rxn_hash->{$bestrxn}->{equation}."\t".$rxn_hash->{$bestrxn}->{name}."\t".$rxn_hash->{$bestrxn}->{gpr}."\t".$rxn_hash->{$bestrxn}->{flux}."\t".$totalflux;
	push(@{$gene_list},$line);
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/Dropbox/workspace/PlantSEED/CCRStudy/".$genome."_gene_flux",$gene_list);