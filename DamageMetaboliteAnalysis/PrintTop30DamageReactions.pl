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
my $infile = "/Users/chenry/workspace/Metabolite repair/Top_30_Reactions.txt";
my $data = Bio::KBase::ObjectAPI::utilities::LOADFILE($infile);
for (my $i=1; $i < @{$data}; $i++) {
	my $array = [split(/\t/,$data->[$i])];
	if (!defined($operators->{$array->[0]})) {
		$operators->{$array->[0]} = {
			name => $array->[0],
			numrxn => 0,
			reactions => {}
		};
	}
	my $rxnarray = [split(/\s+=\>\s+/,$array->[2])];
	my $reactarray = [split(/\s+\+\s+/,$rxnarray->[0])];
	my $prodarray = [split(/\s+\+\s+/,$rxnarray->[1])];
	my $count = @{$reactarray} + @{$prodarray} + 1;
	my $lines = [
	"<tr>",
	"<tr>"
	];
	for (my $j=0; $j < @{$reactarray}; $j++) {
		if ($reactarray->[$j] =~ m/\((\d+)\)\s+(\S+.+)\[(.+)\]/) {
			my $coef = $1;
			my $name = $2;
			my $image = $3;
			$name =~ s/\W//g;
			$lines->[0] .= '<td><img src="'.$image.'" alt="'.$name.'" height="200" width="*"></td>';
			$lines->[1] .= '<td>('.$coef.") ".$name.'</td>';
		}
	}
	$lines->[0] .= "<td>=></td>";
	$lines->[1] .= "<td></td>";
	for (my $j=0; $j < @{$prodarray}; $j++) {
		if ($prodarray->[$j] =~ m/\((\d+)\)\s+(\S+.+)\[(.+)\]/) {
			my $coef = $1;
			my $name = $2;
			my $image = $3;
			$name =~ s/\W//g;
			$lines->[0] .= '<td><img src="'.$image.'" alt="'.$name.'" height="200" width="*"></td>';
			$lines->[1] .= '<td>('.$coef.") ".$name.'</td>';
		}
	}
	$lines->[0] .= "</tr>";
	$lines->[1] .= "</tr>";
	$operators->{$array->[0]}->{numrxn}++;
	$operators->{$array->[0]}->{reactions}->{$array->[1]} = "<table>".join("",@{$lines})."</table>";
}

my $htmlreport = ["<html>","<head>","<script type='text/javascript' src='https://www.google.com/jsapi'></script>","<script type='text/javascript'>","google.load('visualization', '1', {packages:['controls'], callback: drawDashboard});","google.setOnLoadCallback(drawDashboard);"];
push(@{$htmlreport},"function drawDashboard() {var data = new google.visualization.DataTable();");
push(@{$htmlreport},"data.addColumn('string','Operator');");
push(@{$htmlreport},"data.addColumn('string','Reaction count');");
push(@{$htmlreport},"data.addColumn('string','Image');");
push(@{$htmlreport},"data.addRows([");
foreach my $op (keys(%{$operators})) {
	my $row = ['\'<a href="http://bioseed.mcs.anl.gov/~chenry/Website2/'.$op.'.html">'.$op.'</a>\'',"'".$operators->{$op}->{numrxn}."'",'\'<img src="operators/'.$op.'.png" alt="'.$op.'" height="300" width="*">\''];
	push(@{$htmlreport},'['.join(',',@{$row}).'],');
}
push(@{$htmlreport},"]);var filterColumns = [];var tab_columns = [];for (var j = 0, dcols = data.getNumberOfColumns(); j < dcols; j++) {filterColumns.push(j);tab_columns.push(j);}filterColumns.push({type: 'string',calc: function (dt, row) {for (var i = 0, vals = [], cols = dt.getNumberOfColumns(); i < cols; i++) {vals.push(dt.getFormattedValue(row, i));}return vals.join('\\n');}});");
push(@{$htmlreport},"var table = new google.visualization.ChartWrapper({chartType: 'Table',containerId: 'table_div',options: {allowHtml: true,showRowNumber: true,page: 'enable',pageSize: 20},view: {columns: tab_columns}});");
push(@{$htmlreport},"var search_filter = new google.visualization.ControlWrapper({controlType: 'StringFilter',containerId: 'search_div',options: {filterColumnIndex: data.getNumberOfColumns(),matchType: 'any',caseSensitive: false,ui: {label: 'Search data:'}},view: {columns: filterColumns}});");
push(@{$htmlreport},"var dashboard = new google.visualization.Dashboard(document.querySelector('#dashboard_div'));var formatter = new google.visualization.ColorFormat();formatter.addRange(0.5, null, 'red', 'white');");
push(@{$htmlreport},"dashboard.bind([search_filter], [table]);dashboard.draw(data);}</script></head>");
push(@{$htmlreport},"<body><h4>Summary of predicted damage for JCVI minimal genome model</h4><div id='dashboard_div'><table class='columns'><tr><td><div id='search_div'></div></td></tr><tr><td><div id='table_div'></div></td></tr></table></div></body></html>");
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/Website2/Index.html",$htmlreport);

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
	Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/Website2/".$op.".html",$htmlreport);
}
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/CurationTable.txt",$curation_table);
Bio::KBase::ObjectAPI::utilities::PRINTFILE("/Users/chenry/workspace/Metabolite repair/OperatorTable.txt",$operator_table);