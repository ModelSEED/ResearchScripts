#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Spreadsheet::WriteExcel;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $operators = {};
my $infile = "/Users/chenry/workspace/Metabolite repair/iMB155.damage.json";
my $outfile = "/Users/chenry/workspace/Metabolite repair/iMB155.damage.xls";
my $modeldata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($infile)}));
my $cpds = $modeldata->{modelcompounds};
my $cpdhash = {};
chdir "/Users/chenry/workspace/Metabolite repair/compounds/";
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->{id};
	$id =~ s/_[a-z]0//;
	$cpdhash->{$id} = $cpds->[$i];
	if ($id =~ m/cpd\d+/) {
		if (!-e "/Users/chenry/workspace/Metabolite repair/compounds/".$id.".png") {
			system("wget http://minedatabase.mcs.anl.gov/compound_images/ModelSEED/".$id.".png");
		}
	}
}
exit;
my $rxns = $modeldata->{modelreactions};
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	if (!defined($operators->{$rxn->{reference}})) {
		$operators->{$rxn->{reference}} = {
			name => $rxn->{reference},
			numcpd => 0,
			numrxn => 0,
			reactants => {},
			products => {},
			primary_reactants => {},
			primary_products => {},
			cofactor_reactants => {},
			cofactor_products => {},
			reactions => {}
		};
	}
	$operators->{$rxn->{reference}}->{numrxn}++;
	for (my $j=0; $j < @{$rxn->{modelReactionReagents}}; $j++) {
		my $id = $rxn->{modelReactionReagents}->[$j]->{modelcompound_ref};
		$id =~ s/.+\///g;
		if ($rxn->{modelReactionReagents}->[$j]->{modelcompound_ref} =~ m/(cpd\d+)/) {
			$id = $1;
		} elsif ($rxn->{modelReactionReagents}->[$j]->{modelcompound_ref} =~ m/(pkc\d+)/) {
			$id = $1;
		}
		if ($rxn->{modelReactionReagents}->[$j]->{coefficient} < 0) {
			if (!defined($operators->{$rxn->{reference}}->{reactants}->{$id})) {
				$operators->{$rxn->{reference}}->{reactants}->{$id} = 0;
			}
			$operators->{$rxn->{reference}}->{reactants}->{$id}++;
		} else {
			if (!defined($operators->{$rxn->{reference}}->{products}->{$id})) {
				$operators->{$rxn->{reference}}->{products}->{$id} = 0;
			}
			$operators->{$rxn->{reference}}->{products}->{$id}++;
		}
	}
}
my $maxcofactor = 0;
my $maxprimary = 0;
foreach my $op (keys(%{$operators})) {
	foreach my $cpd (keys(%{$operators->{$op}->{reactants}})) {
		my $fraction = $operators->{$op}->{reactants}->{$cpd}/$operators->{$op}->{numrxn};
		if ($fraction < 0.5 || $cpd !~ m/cpd/) {
			if (!defined($operators->{$op}->{primary_reactants}->{$cpd})) {
				$operators->{$op}->{primary_reactants}->{$cpd} = 0;
			}
			$operators->{$op}->{primary_reactants}->{$cpd}++;
		} else {
			if (!defined($operators->{$op}->{cofactor_reactants}->{$cpd})) {
				$operators->{$op}->{cofactor_reactants}->{$cpd} = 0;
			}
			$operators->{$op}->{cofactor_reactants}->{$cpd}++;
		}
	}
	foreach my $cpd (keys(%{$operators->{$op}->{products}})) {
		my $fraction = $operators->{$op}->{products}->{$cpd}/$operators->{$op}->{numrxn};
		if ($fraction < 0.5 || $cpd !~ m/cpd/) {
			if (!defined($operators->{$op}->{primary_products}->{$cpd})) {
				$operators->{$op}->{primary_products}->{$cpd} = 0;
			}
			$operators->{$op}->{primary_products}->{$cpd}++;
		} else {
			if (!defined($operators->{$op}->{cofactor_products}->{$cpd})) {
				$operators->{$op}->{cofactor_products}->{$cpd} = 0;
			}
			$operators->{$op}->{cofactor_products}->{$cpd}++;
		}
	}
	$operators->{$op}->{numcpd} = keys(%{$operators->{$op}->{primary_reactants}});
	my $numcofprod = keys(%{$operators->{$op}->{cofactor_products}});
	my $numcofreact = keys(%{$operators->{$op}->{cofactor_reactants}});
	if (($numcofprod + $numcofreact) > $maxcofactor) {
		$maxcofactor = ($numcofprod + $numcofreact);
	}
	$numcofprod = keys(%{$operators->{$op}->{primary_reactants}});
	$numcofreact = keys(%{$operators->{$op}->{primary_products}});
	if (($numcofprod + $numcofreact) > $maxprimary) {
		$maxprimary = ($numcofprod + $numcofreact);
	}
}

my $wkbk = Spreadsheet::WriteExcel->new($outfile);
my $sheet = $wkbk->add_worksheet("Overview");
my $headings = ["Operator","Compound count","Reaction count"];
for (my $i=0; $i < $maxcofactor; $i++) {
	my $index = $i+1;
	push(@{$headings},"Cofactor ".$index);
}
push(@{$headings},"Image");
$sheet->write_row(0,0,$headings);
my $row = 1;
chdir "/Users/chenry/workspace/Metabolite repair/operators/";
foreach my $op (keys(%{$operators})) {
	if (!-e "/Users/chenry/workspace/Metabolite repair/operators/".$op.".png") {
		system("wget http://webfba.chem-eng.northwestern.edu/MINE_imgs/op_images/".$op.".png");
	}
	my $newrow = [$op,$operators->{$op}->{numcpd},$operators->{$op}->{numrxn}];
	for (my $i=0; $i < $maxcofactor; $i++) {
		push(@{$newrow},"");
	}
	my $cofcount = 0;
	foreach my $cpd (keys(%{$operators->{$op}->{cofactor_reactants}})) {
		$newrow->[$cofcount+3] = "(R) ".$cpdhash->{$cpd}->{name};
		$cofcount++;
	}
	foreach my $cpd (keys(%{$operators->{$op}->{cofactor_products}})) {
		$newrow->[$cofcount+3] = "(P) ".$cpdhash->{$cpd}->{name};
		$cofcount++;
	}
	$sheet->write_row($row,0,$newrow);
	my $column = "K";
	$sheet->insert_image($column.$row,"/Users/chenry/workspace/Metabolite repair/Operators/".$op.".png");
	$row++;
}