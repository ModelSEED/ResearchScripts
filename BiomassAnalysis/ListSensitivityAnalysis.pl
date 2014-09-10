#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $workspace = $ARGV[0];
my $output = get_ws_objects_list($workspace,"KBaseFBA.ReactionSensitivityAnalysis");
my $tbl = [];
my $opt = {};
for (my $i=0; $i < @{$output};$i++) {
    my $r = $output->[$i];
    my $size = $r->[9]+0;
    if (defined($opt->{megabytes})) {
	$size = int(($size/1048576)*1000+0.5)/1000; # convert to MB, rounded to three decimals
    }
    push(@{$tbl},[$r->[0],$r->[1],$r->[4],$r->[2],$r->[6],$r->[7],$r->[5],$r->[3],$size]);
}
my $sizeHeader = 'Size(bytes)';
if (defined($opt->{megabytes})) {
	$sizeHeader = 'Size(MB)';
}
my $table = Text::Table->new(
	'ID', 'ObjName', 'Vers', 'Type','WSID','WS','Last_modby','Moddate',$sizeHeader
	);
my @sorted_tbl = @$tbl;
if (defined($opt->{column})) {
	if ($opt->{column}==9) {
		#size is numeric, so sort numerically, largest first
		@sorted_tbl = sort { $b->[$opt->{column}-1] <=> $a->[$opt->{column}-1] } @sorted_tbl;
	} elsif ( $opt->{column}==1 || $opt->{column}==3 || $opt->{column}==5) {
		#ids and version numbers are numeric, so sort numerically, largest last
		@sorted_tbl = sort { $a->[$opt->{column}-1] <=> $b->[$opt->{column}-1] } @sorted_tbl;
	} else {
		@sorted_tbl = sort { $a->[$opt->{column}-1] cmp $b->[$opt->{column}-1] } @sorted_tbl;
	}
}
# splice out the first n if limit is set
if (defined($opt->{limit})) {
	@sorted_tbl=splice(@sorted_tbl,0,$opt->{limit});
}

$table->load(@sorted_tbl);
print $table;
	