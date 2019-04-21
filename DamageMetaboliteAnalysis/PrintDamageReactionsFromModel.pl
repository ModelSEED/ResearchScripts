#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Spreadsheet::WriteExcel;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $infile = "/Users/chenry/Dropbox/workspace/Metabolite repair/iMB155.damage.json";
my $modeldata = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE($infile)}));
my $cpds = $modeldata->{modelcompounds};
my $cpdhash = {};
for (my $i=0; $i < @{$cpds}; $i++) {
	my $id = $cpds->[$i]->{id};
	$id =~ s/_c0//;
	$cpdhash->{$id} = $cpds->[$i];
	if (!defined($cpdhash->{$id}->{name})) {
		$cpdhash->{$id}->{name} = $id;
	}
	$cpdhash->{$id}->{name} =~ s/_c0//;
}
my $operators = {};
my $rxns = $modeldata->{modelreactions};
for (my $i=0; $i < @{$rxns}; $i++) {
	my $rxn = $rxns->[$i];
	$rxn->{id} =~ s/_c0//;
	if (!defined($operators->{$rxn->{reference}})) {
		$operators->{$rxn->{reference}} = {
			name => $rxn->{reference},
			numrxn => 0,
			reactions => {}
		};
	}
	$operators->{$rxn->{reference}}->{numrxn}++;
	my $lines = [
	"<tr>",
	"<tr>"
	];
	my $prodlines = ["",""];
	for (my $j=0; $j < @{$rxn->{modelReactionReagents}}; $j++) {
		my $id = $rxn->{modelReactionReagents}->[$j]->{modelcompound_ref};
		$id =~ s/.+\///g;
		if ($rxn->{modelReactionReagents}->[$j]->{modelcompound_ref} =~ m/(cpd\d+)/) {
			$id = $1;
		} elsif ($rxn->{modelReactionReagents}->[$j]->{modelcompound_ref} =~ m/(pkc\d+)/) {
			$id = $1;
		}
		if ($rxn->{modelReactionReagents}->[$j]->{coefficient} < 0) {
			$lines->[0] .= '<td><img src="compounds/'.$id.'.png" alt="'.$cpdhash->{$id}->{name}.'" height="200" width="*"></td>';
			$lines->[1] .= '<td>('.-1*$rxn->{modelReactionReagents}->[$j]->{coefficient}.") ".$cpdhash->{$id}->{name}.'</td>';
		} else {
			$prodlines->[0] .= '<td><img src="compounds/'.$id.'.png" alt="'.$cpdhash->{$id}->{name}.'" height="200" width="*"></td>';
			$prodlines->[1] .= '<td>('.$rxn->{modelReactionReagents}->[$j]->{coefficient}.") ".$cpdhash->{$id}->{name}.'</td>';
		}
	}
	$lines->[0] .= "<td>=></td>".$prodlines->[0]."</tr>";
	$lines->[1] .= "<td></td>".$prodlines->[1]."</tr>";
	$operators->{$rxn->{reference}}->{reactions}->{$rxn->{id}} = "<table>".join("",@{$lines})."</table>";
}	
my $htmlreport = ["<html>","<head>","<script type='text/javascript' src='https://www.google.com/jsapi'></script>","<script type='text/javascript'>","google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});","google.setOnLoadCallback(drawDashboard);"];
push(@{$htmlreport},"function drawDashboard() {var data = new google.visualization.DataTable();");
push(@{$htmlreport},"data.addColumn('string','Operator');");
push(@{$htmlreport},"data.addColumn('string','Reaction count');");
push(@{$htmlreport},"data.addColumn('string','Image');");
push(@{$htmlreport},"data.addRows([");
foreach my $op (keys(%{$operators})) {
	my $row = ['\'<a href="http://bioseed.mcs.anl.gov/~chenry/Website/'.$op.'.html">'.$op.'</a>\'',"'".$operators->{$op}->{numrxn}."'",'\'<img src="operators/'.$op.'.png" alt="'.$op.'" height="300" width="*">\''];
	push(@{$htmlreport},'['.join(',',@{$row}).'],');
}
push(@{$htmlreport},"]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});");
push(@{$htmlreport},"var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});");
push(@{$htmlreport},"var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});");
push(@{$htmlreport},"var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');");
push(@{$htmlreport},"dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>");
push(@{$htmlreport},"<body><h4>Summary of predicted damage for JCVI minimal genome model</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>");
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/Website/Index.html",$htmlreport);

my $curation_table = [];
my $operator_table = [];
foreach my $op (keys(%{$operators})) {
	push(@{$operator_table},$op);
	$htmlreport = ["<html>","<head>","<script type='text/javascript' src='https://www.google.com/jsapi'></script>","<script type='text/javascript'>","google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});","google.setOnLoadCallback(drawDashboard);"];
	push(@{$htmlreport},"function drawDashboard() {var data = new google.visualization.DataTable();");
	push(@{$htmlreport},"data.addColumn('string','Reaction');");
	push(@{$htmlreport},"data.addColumn('string','Equation');");
	push(@{$htmlreport},"data.addRows([");
	my $index = 1;
	foreach my $rxn (keys(%{$operators->{$op}->{reactions}})) {
		push(@{$curation_table},$op." ".$index."\t".$rxn);
		my $row = ["'".$rxn."'","'".$operators->{$op}->{reactions}->{$rxn}."'"];
		push(@{$htmlreport},'['.join(',',@{$row}).'],');
		$index++;
	}
	push(@{$htmlreport},"]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});");
	push(@{$htmlreport},"var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});");
	push(@{$htmlreport},"var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});");
	push(@{$htmlreport},"var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');");
	push(@{$htmlreport},"dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>");
	push(@{$htmlreport},"<body><h4>Reactions for operator ".$op."</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>");
	Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/Website/".$op.".html",$htmlreport);
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/CurationTable.txt",$curation_table);
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/OperatorTable.txt",$operator_table);