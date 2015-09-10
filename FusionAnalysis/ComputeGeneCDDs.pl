#!/usr/bin/perl

use strict;
use warnings;
use ModelSEED::Client::SAP;
use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $directory = $ARGV[0];

open(my $fhh, "<", $directory."CDD-Data.txt");
my $line = <$fhh>;
my $cdddata = {};
while ($line = <$fhh>) {
	my $array = [split(/\t/,$line)];
	$cdddata->{$array->[0]} = {
		id => $array->[0],
		len => $array->[1],
		name => $array->[2],
		genes => $array->[3],
		singlegenes => $array->[4],
		longgenes => $array->[5]
	};
}
close($fhh);

my $genes = [qw(
fig|226186.1.peg.1442
fig|525146.3.peg.1509  
fig|555079.3.peg.2426
fig|233413.1.peg.32 
fig|498211.3.peg.425  
fig|267747.1.peg.1824
fig|528347.5.peg.721  
fig|75379.4.peg.581
fig|469382.4.peg.3193  
fig|4932.3.peg.877
fig|272557.1.peg.1396  
fig|313606.3.peg.7087
fig|563038.3.peg.1439  
fig|592015.5.peg.1366  
fig|391774.5.peg.1943  
fig|373903.5.peg.523  
fig|1001739.3.peg.2350
fig|188626.3.peg.1493 
fig|247156.1.peg.747
fig|59931.3.peg.1325  
fig|298653.4.peg.1987
fig|639282.3.peg.2135
fig|591167.6.peg.6429
fig|546274.4.peg.1959
fig|471853.5.peg.728
fig|525246.3.peg.1726  
fig|180281.4.peg.761  
fig|546271.3.peg.778  
fig|439235.3.peg.5130  
fig|411465.10.peg.1479  
fig|405948.11.peg.5492 
fig|1127122.3.peg.1544
fig|10090.3.peg.25319  
fig|3702.11.peg.9804
fig|4896.1.peg.4275  
fig|376686.6.peg.1713  
fig|500635.8.peg.2116
fig|309798.3.peg.1142 
fig|272561.1.peg.628 
fig|521095.6.peg.109
fig|381764.6.peg.1714 
fig|382464.3.peg.4351
fig|431947.6.peg.1037
fig|251221.1.peg.2647
fig|391603.3.peg.524
fig|233150.3.peg.42
fig|478801.5.peg.881
fig|272634.1.peg.300
fig|471853.5.peg.717
fig|330779.3.peg.1576
fig|349741.3.peg.1549  
fig|264201.15.peg.1181
fig|203267.1.peg.586
fig|395965.4.peg.2584 
fig|481448.7.peg.1889 
fig|290400.10.peg.2898 
fig|765698.3.peg.4766
fig|378806.7.peg.2871  
fig|367737.4.peg.1188
fig|519442.4.peg.1288
fig|390333.6.peg.218
fig|469381.4.peg.1963
fig|177437.4.peg.854 
fig|257309.1.peg.1719 
fig|138119.3.peg.3934
fig|3702.11.peg.16955
fig|4932.3.peg.5208
fig|668336.4.peg.1151
fig|500633.7.peg.413
fig|497965.6.peg.4479Ê 
fig|161528.3.peg.1164
fig|314265.3.peg.941
fig|63737.4.peg.3651
fig|438753.3.peg.3255  
fig|272559.3.peg.965
fig|357808.3.peg.923  
fig|101510.15.peg.5449
fig|292805.3.peg.329 
fig|321955.3.peg.1718
fig|223926.6.peg.2972
fig|452863.6.peg.2446 
fig|69014.3.peg.309  
fig|100226.1.peg.4062
fig|426368.9.peg.603  
fig|240015.3.peg.2233
fig|479435.6.peg.314  
fig|644966.3.peg.1460 
fig|83333.1.peg.25
fig|608534.3.peg.920  
fig|708616.3.peg.127  
fig|64091.1.peg.1812  
fig|316067.3.peg.130  
fig|330779.3.peg.297
fig|4932.3.peg.5667  
fig|3702.11.peg.14378
fig|3702.11.peg.466
fig|3702.11.peg.8629 
fig|572480.3.peg.1028
fig|452638.3.peg.1662
fig|64091.1.peg.2192  
fig|395495.3.peg.2329  
fig|590998.5.peg.1658
fig|273057.1.peg.2
fig|411474.6.peg.1072  
fig|866775.3.peg.61 
fig|243277.1.peg.62
fig|264198.3.peg.721
fig|234826.3.peg.438
fig|196164.1.peg.1591
fig|515620.4.peg.2275  
fig|96561.3.peg.1759 
fig|717606.6.peg.1174
fig|187303.17.peg.2469
fig|469616.3.peg.1833  
fig|326426.4.peg.1408 
fig|342108.5.peg.2398  
fig|338963.3.peg.1505
fig|479436.6.peg.307  
fig|471875.6.peg.815
fig|431947.6.peg.153
fig|471854.4.peg.5604
fig|521095.6.peg.839
fig|103690.1.peg.3826
fig|554065.3.peg.6386
fig|4896.1.peg.3261 
fig|3702.11.peg.2993
fig|269800.4.peg.2728
fig|178306.1.peg.2190  
fig|272635.1.peg.335
fig|272568.11.peg.1069
fig|178306.1.peg.2192  
fig|1148.1.peg.2723
fig|523794.5.peg.1959 
fig|83333.1.peg.4298
fig|243232.1.peg.949
fig|1148.1.peg.1732
fig|585501.3.peg.1201
fig|362976.10.peg.504  
fig|240016.6.peg.2924
fig|391625.5.peg.3242
fig|3702.11.peg.4621
fig|3702.11.peg.22518
fig|3702.11.peg.4649
fig|7227.3.peg.13624
fig|9606.3.peg.25573
fig|266779.9.peg.167
fig|219305.4.peg.3726
fig|177439.1.peg.1605  
fig|632772.3.peg.4758  
fig|3702.7.peg.25891	
)];

my $sapsvr = ModelSEED::Client::SAP->new();
my $functions = $sapsvr->ids_to_functions({-ids => $genes});
my $sequences = $sapsvr->fids_to_proteins({-ids => $genes,-sequence => 1});
my $output = "";
foreach my $gene (keys(%{$sequences})) {
	$output .= ">".$gene."\n".$sequences->{$gene}."\n";
}
open(my $fh, ">", $directory."genes.fasta");
print $fh $output;
close ($fh);

system("rpsblast -query ".$directory."genes.fasta -evalue 0.01 -seg no -outfmt 6 -db /disks/olson/cdd/data/Cdd -out ".$directory."genes.cdd");

open(my $fhhhh, ">", $directory."genes.cdd.out");
open(my $fhhh, "<", $directory."genes.cdd");
print $fhhhh "Gene	Length	CDD	Start	Stop	Identity	Function	Alignlength	CDD name	Protein	CDDStart	CDDStop	E value	SEEDID\n";
while (my $line = <$fhhh>) {
	my $array = [split(/\t/,$line)];
	my $itemarray = [split(/\|/,$array->[1])];
	print $fhhhh $array->[0]."\t".3*length($sequences->{$array->[0]})."\t".$itemarray->[2]."\t".
		$array->[6]."\t".$array->[7]."\t".$array->[2]."\t".$functions->{$array->[0]}."\t".
		$array->[3]."\t".$cdddata->{$itemarray->[2]}."\t".Digest::MD5::md5_hex($sequences->{$array->[0]})."\t".$array->[8]."\t".$array->[9]."\t".
		$array->[10]."\t".$genes->{$array->[0]}->{id}."\n";
}
close($fhhh);
close($fhhhh);