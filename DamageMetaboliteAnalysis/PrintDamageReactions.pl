#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Spreadsheet::WriteExcel;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $letters = [qw(
A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
)];
my $operators = {};
my $infile = "/Users/christopherhenry/Dropbox/workspace/Metabolite repair/iMB155.damage.json";
my $outfile = "/Users/christopherhenry/Dropbox/workspace/Metabolite repair/iMB155.damage.xls";
my $modeldata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($infile)}));
my $cpds = $modeldata->{modelcompounds};
my $cpdhash = {};
chdir "/Users/christopherhenry/Dropbox/workspace/Metabolite repair/compounds/";
for (my $i=0; $i < @{$cpds}; $i++) {
	$cpds->[$i]->{id} =~ s/_[a-z]0//;
	$cpds->[$i]->{name} =~ s/_[a-z]0//;
	#print $cpds->[$i]->{id}."\t".$cpds->[$i]->{smiles}."\n";
	$cpdhash->{$cpds->[$i]->{id}} = $cpds->[$i];
	$cpdhash->{$cpds->[$i]->{id}}->{count} = 0;
	if ($cpds->[$i]->{id} =~ m/cpd\d+/) {
		if (!-e "/Users/christopherhenry/Dropbox/workspace/Metabolite repair/compounds/".$cpds->[$i]->{id}.".png") {
			#system("wget http://minedatabase.mcs.anl.gov/compound_images/ModelSEED/".$cpds->[$i]->{id}.".png");
		}
	}
}
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
		$cpdhash->{$id}->{count}++;
		if ($rxn->{modelReactionReagents}->[$j]->{coefficient} < 0) {
			if (!defined($operators->{$rxn->{reference}}->{reactants}->{$id})) {
				$operators->{$rxn->{reference}}->{reactants}->{$id}  = 0;
			}
			$operators->{$rxn->{reference}}->{reactants}->{$id}++;
			$operators->{$rxn->{reference}}->{reactions}->{$rxn->{id}}->{reactants}->{$id} = 1;
		} else {
			if (!defined($operators->{$rxn->{reference}}->{products}->{$id})) {
				$operators->{$rxn->{reference}}->{products}->{$id}  = 0;
			}
			$operators->{$rxn->{reference}}->{products}->{$id}++;
			$operators->{$rxn->{reference}}->{reactions}->{$rxn->{id}}->{products}->{$id} = 1;
		}
	}
}
foreach my $op (keys(%{$operators})) {
	foreach my $cpd (keys(%{$operators->{$op}->{reactants}})) {
		my $fraction = $operators->{$op}->{reactants}->{$cpd}/$operators->{$op}->{numrxn};
		if ($fraction < 1 || $cpd !~ m/cpd/ || ($operators->{$op}->{numrxn} == 1 && $cpdhash->{$cpd}->{count} < 5)) {
			$cpdhash->{$cpd}->{ops}->{$op} = "R";
			if (!defined($operators->{$op}->{primary_reactants}->{$cpd})) {
				$operators->{$op}->{primary_reactants}->{$cpd} = 0;
			}
			$operators->{$op}->{primary_reactants}->{$cpd}++;
		} else {
			$cpdhash->{$cpd}->{ops}->{$op} = "C";
			if (!defined($operators->{$op}->{cofactor_reactants}->{$cpd})) {
				$operators->{$op}->{cofactor_reactants}->{$cpd} = 0;
			}
			$operators->{$op}->{cofactor_reactants}->{$cpd}++;
		}
	}
	foreach my $cpd (keys(%{$operators->{$op}->{products}})) {
		my $fraction = $operators->{$op}->{products}->{$cpd}/$operators->{$op}->{numrxn};
		if ($fraction < 1 || $cpd !~ m/cpd/ || ($operators->{$op}->{numrxn} == 1 && $cpdhash->{$cpd}->{count} < 5)) {
			$cpdhash->{$cpd}->{ops}->{$op} = "R";
			if (!defined($operators->{$op}->{primary_products}->{$cpd})) {
				$operators->{$op}->{primary_products}->{$cpd} = 0;
			}
			$operators->{$op}->{primary_products}->{$cpd}++;
		} else {
			$cpdhash->{$cpd}->{ops}->{$op} = "C";
			if (!defined($operators->{$op}->{cofactor_products}->{$cpd})) {
				$operators->{$op}->{cofactor_products}->{$cpd} = 0;
			}
			$operators->{$op}->{cofactor_products}->{$cpd}++;
		}
	}
	$operators->{$op}->{numcpd} = keys(%{$operators->{$op}->{primary_reactants}});
}

my $htmlreport = ["<html>","<head>","<script type='text/javascript' src='https://www.google.com/jsapi'></script>","<script type='text/javascript'>","google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});","google.setOnLoadCallback(drawDashboard);"];
push(@{$htmlreport},"function drawDashboard() {var data = new google.visualization.DataTable();");
push(@{$htmlreport},"data.addColumn('string','Operator');");
push(@{$htmlreport},"data.addColumn('string','Compound count');");
push(@{$htmlreport},"data.addColumn('string','Reaction count');");
push(@{$htmlreport},"data.addColumn('string','Cofactors');");
push(@{$htmlreport},"data.addColumn('string','Image');");
push(@{$htmlreport},"data.addRows([");
foreach my $op (keys(%{$operators})) {
	my $row = ['\'<a href="http://bioseed.mcs.anl.gov/~chenry/Website/'.$op.'.html">'.$op.'</a>\'',"'".$operators->{$op}->{numcpd}."'","'".$operators->{$op}->{numrxn}."'","",'\'<img src="operators/'.$op.'.png" alt="'.$op.'" height="300" width="*">\''];
	$row->[3] = "'";
	foreach my $cpd (keys(%{$operators->{$op}->{cofactor_reactants}})) {
		if (length($row->[3]) > 1) {
			$row->[3] .= "<br>";
		}
		$row->[3] .= "(R) ".$cpdhash->{$cpd}->{name};
	}
	foreach my $cpd (keys(%{$operators->{$op}->{cofactor_products}})) {
		if (length($row->[3]) > 1) {
			$row->[3] .= "<br>";
		}
		$row->[3] .= "(P) ".$cpdhash->{$cpd}->{name};
	}
	$row->[3] .= "'";
	push(@{$htmlreport},'['.join(',',@{$row}).'],');
}
push(@{$htmlreport},"]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});");
push(@{$htmlreport},"var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});");
push(@{$htmlreport},"var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});");
push(@{$htmlreport},"var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');");
push(@{$htmlreport},"dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>");
push(@{$htmlreport},"<body><h4>Summary of predicted damage for JCVI minimal genome model</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>");
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/christopherhenry/Dropbox/workspace/Metabolite repair/Website/Index.html",$htmlreport);

foreach my $op (keys(%{$operators})) {
	$htmlreport = ["<html>","<head>","<script type='text/javascript' src='https://www.google.com/jsapi'></script>","<script type='text/javascript'>","google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});","google.setOnLoadCallback(drawDashboard);"];
	push(@{$htmlreport},"function drawDashboard() {var data = new google.visualization.DataTable();");
	push(@{$htmlreport},"data.addColumn('string','Reaction');");
	push(@{$htmlreport},"data.addColumn('string','Primary reactant');");
	push(@{$htmlreport},"data.addColumn('string','Cofactors');");
	push(@{$htmlreport},"data.addColumn('string','Primary product');");
	push(@{$htmlreport},"data.addRows([");
	foreach my $rxn (keys(%{$operators->{$op}->{reactions}})) {
		my $row = ["'".$rxn."'","'","'","'"];
		foreach my $cpd (keys(%{$operators->{$op}->{reactions}->{$rxn}->{reactants}})) {
			if ($cpdhash->{$cpd}->{ops}->{$op} eq "R") {
				if (length($row->[1]) > 1) {
					$row->[1] .= "<br>";
				}
				$row->[1] .= '<img src="compounds/'.$cpd.'.png" alt="'.$op.'" height="200" width="*"><br>'.$cpdhash->{$cpd}->{name};
			} else {
				if (length($row->[2]) > 1) {
					$row->[2] .= "<br>";
				}
				$row->[2] .= "(R) ".$cpdhash->{$cpd}->{name};
			}
		}
		foreach my $cpd (keys(%{$operators->{$op}->{reactions}->{$rxn}->{products}})) {
			if ($cpdhash->{$cpd}->{ops}->{$op} eq "R") {
				if (length($row->[3]) > 1) {
					$row->[3] .= "<br>";
				}
				$row->[3] .= '<img src="compounds/'.$cpd.'.png" alt="'.$op.'" height="200" width="*"><br>'.$cpdhash->{$cpd}->{name};
			} else {
				if (length($row->[2]) > 1) {
					$row->[2] .= "<br>";
				}
				$row->[2] .= "(P) ".$cpdhash->{$cpd}->{name};
			}
		}
		$row->[1] .= "'";
		$row->[2] .= "'";
		$row->[3] .= "'";
		push(@{$htmlreport},'['.join(',',@{$row}).'],');
	}
	push(@{$htmlreport},"]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});");
	push(@{$htmlreport},"var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});");
	push(@{$htmlreport},"var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});");
	push(@{$htmlreport},"var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');");
	push(@{$htmlreport},"dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>");
	push(@{$htmlreport},"<body><h4>Reactions for operator ".$op."</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>");
	Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/christopherhenry/Dropbox/workspace/Metabolite repair/Website/".$op.".html",$htmlreport);
}
	