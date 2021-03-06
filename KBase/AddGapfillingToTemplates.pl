#!/usr/bin/perl
   
use strict;
use Data::Dumper;
use JSON::XS;
use DateTime;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(fbaws get_fba_client runFBACommand universalFBAScriptCode );
use Bio::KBase::ObjectAPI::KBaseStore;

$|=1;

my $store = Bio::KBase::ObjectAPI::KBaseStore->new({
	provenance => [],
	workspace => get_ws_client()
});

my $unbalrxn = [qw(
rxn04770_c0_<
rxn02777_c0_>
rxn00750_c0_<
rxn06079_c0_<
rxn06219_c0_>
rxn06629_c0_>
rxn06850_c0_>
rxn07172_c0_>
rxn07180_c0_>
rxn07351_c0_>
rxn07581_c0_>
rxn07649_c0_>
rxn10083_c0_<
rxn10430_c0_>
rxn10446_c0_>
rxn10458_c0_>
rxn10459_c0_>
rxn10460_c0_>
rxn10461_c0_>
rxn10486_c0_>
rxn10503_c0_>
rxn10532_c0_>
rxn10546_c0_>
rxn10549_c0_>
rxn10566_c0_>
rxn10582_c0_>
rxn10642_c0_<
rxn10652_c0_<
rxn10659_c0_<
rxn10685_c0_>
rxn10686_c0_>
rxn10687_c0_>
rxn10697_c0_>
rxn10698_c0_>
rxn10710_c0_>
rxn10724_c0_>
rxn10726_c0_>
rxn10736_c0_>
rxn10758_c0_>
rxn10759_c0_>
rxn10765_c0_>
rxn10773_c0_>
rxn10780_c0_>
rxn10808_c0_<
rxn11197_c0_<
rxn12899_c0_>
rxn13472_c0_<
rxn13505_c0_>
rxn13515_c0_>
rxn13517_c0_<
rxn13674_c0_>
rxn13675_c0_>
)];

my $gnegtemp = $store->get_object("KBaseTemplateModels/GramNegModelTemplate");
my $biocpds = $gnegtemp->templateBiomasses()->[0]->templateBiomassComponents();
my $cpdhash = {};
foreach my $biocpd (@{$biocpds}) {
	$cpdhash->{$biocpd->compound()->id()} = $biocpd;	
}
my $temprxns = $gnegtemp->templateReactions();
my $rxnhash = {};
foreach my $rxn (@{$temprxns}) {
	$rxnhash->{$rxn->reaction()->id()} = 1;	
}
my $count = @{$temprxns};
for (my $i=0; $i < @{$unbalrxn}; $i++) {
	if ($unbalrxn->[$i] =~ m/(rxn\d+)_c0_(.)/) {
		my $rxn = $1;
		my $dir = $2;
		if (!defined($rxnhash->{$rxn})) {
			$gnegtemp->add("templateReactions",{
				id => $gnegtemp->id().".temprxn.".$count,
		    	reaction_ref => "489/6/1/reactions/id/".$rxn,
		    	compartment_ref => "489/6/1/compartments/id/c",
		    	complex_refs => [],
		    	direction => $dir,
		    	type => "universal"
			});
			$count++;
		}
	}
}
$count = @{$biocpds};
my $gpostemp = $store->get_object("KBaseTemplateModels/GramPosModelTemplate");
$biocpds = $gpostemp->templateBiomasses()->[0]->templateBiomassComponents();
foreach my $biocpd (@{$biocpds}) {
	if (!defined($cpdhash->{$biocpd->compound()->id()})) {
		$count++;
		$biocpd->id($gnegtemp->templateBiomasses()->[0]->id().".tempbiocomp.".$count);
		$gnegtemp->templateBiomasses()->[0]->add("templateBiomassComponents",$biocpd);
	}
}
my $compounds = [qw(
cpd00004
cpd00005
cpd00013
cpd00018
cpd00020
cpd00022
cpd00024
cpd00026
cpd00027
cpd00031
cpd00036
cpd00037
cpd00040
cpd00043
cpd00046
cpd00050
cpd00061
cpd00070
cpd00078
cpd00079
cpd00083
cpd00086
cpd00089
cpd00091
cpd00096
cpd00097
cpd00100
cpd00101
cpd00102
cpd00104
cpd00105
cpd00106
cpd00108
cpd00117
cpd00122
cpd00125
cpd00126
cpd00130
cpd00134
cpd00137
cpd00138
cpd00151
cpd00155
cpd00164
cpd00169
cpd00182
cpd00184
cpd00186
cpd00202
cpd00206
cpd00214
cpd00216
cpd00235
cpd00244
cpd00249
cpd00274
cpd00276
cpd00279
cpd00284
cpd00285
cpd00286
cpd00288
cpd00293
cpd00294
cpd00296
cpd00298
cpd00304
cpd00305
cpd00311
cpd00327
cpd00330
cpd00396
cpd00472
cpd00482
cpd00485
cpd00504
cpd00516
cpd00519
cpd00523
cpd00541
cpd00644
cpd00731
cpd00755
cpd00773
cpd00794
cpd00832
cpd00834
cpd00842
cpd00875
cpd00895
cpd00982
cpd01080
cpd01122
cpd01270
cpd01476
cpd01741
cpd02113
cpd02197
cpd02246
cpd02333
cpd02557
cpd02574
cpd02611
cpd02817
cpd02993
cpd03217
cpd03425
cpd03444
cpd03445
cpd03453
cpd03487
cpd03491
cpd03495
cpd03517
cpd03671
cpd03802
cpd03847
cpd03848
cpd03850
cpd08995
cpd09680
cpd11312
cpd11313
cpd11430
cpd11431
cpd11433
cpd11436
cpd11438
cpd11440
cpd11466
cpd11468
cpd11474
cpd11476
cpd11652
cpd11715
cpd11746
cpd11825
cpd12371
cpd12813
cpd12816
cpd12836
cpd14231
cpd14612
cpd14630
cpd15237
cpd15239
cpd15240
cpd15241
cpd15268
cpd15269
cpd15294
cpd15297
cpd15353
cpd15378
cpd15425
cpd15426
cpd15428
cpd15429
cpd15430
cpd15431
cpd15433
cpd15434
cpd15435
cpd15436
cpd15486
cpd15499
cpd15501
cpd15503
cpd15505
cpd15506
cpd15508
cpd15511
cpd15514
cpd15521
cpd15522
cpd15523
cpd15524
cpd15525
cpd15526
cpd15527
cpd15528
cpd15529
cpd15531
cpd15532
cpd15534
cpd15535
cpd15536
cpd15537
cpd15538
cpd15539
cpd15541
cpd15545
cpd15555
cpd15556
cpd15557
cpd15561
cpd15601
cpd15607
cpd15609
cpd15647
cpd15648
cpd15697
cpd15698
cpd15699
cpd15700
cpd15724
cpd15725
cpd15726
cpd15727
cpd15782
cpd15783
cpd15784
cpd15785
cpd15786
cpd15787
cpd15788
cpd15789
cpd15790
cpd15817
cpd15818
cpd15819
cpd15820
cpd15834
cpd15858
cpd15859
cpd15860
cpd15861
cpd15862
cpd15863
cpd15866
cpd15868
cpd15869
cpd15870
cpd15871
cpd15872
cpd15880
cpd15884
cpd15915
cpd15916
cpd15917
cpd15918
cpd15919
cpd15920
cpd15921
cpd15922
cpd15923
cpd15924
cpd15925
cpd15929
cpd15939
cpd15952
cpd15968
cpd15971
cpd15977
cpd15982
cpd15990
cpd15992
cpd16001
cpd16004
cpd16006
cpd16007
cpd16009
cpd16010
cpd16013
cpd16022
cpd16035
cpd16040
cpd16044
cpd16045
cpd16046
cpd16047
cpd16051
cpd16053
cpd16969
cpd17018
cpd17021
cpd17029
)];
foreach my $biocpd (@{$compounds}) {
	if (!defined($cpdhash->{$biocpd})) {
		$count++;
		$gnegtemp->templateBiomasses()->[0]->add("templateBiomassComponents",{
			id => $gnegtemp->templateBiomasses()->[0]->id().".tempbiocomp.".$count,
			class => "cofactor",
			compound_ref => "489/6/1/compounds/id/".$biocpd,
			compartment_ref => "489/6/1/compartments/id/c",
			coefficientType => "MOLSPLIT",
			coefficient => -1,
			linked_compound_refs => [],
			link_coefficients => []
		});
	}
}
$gnegtemp->save("chenrydemo/FullBiomassTemplate");
1;