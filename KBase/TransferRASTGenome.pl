#!/usr/bin/perl -w

use strict;
use Bio::ModelSEED::MSSeedSupportServer::MSSeedSupportClient;
use Bio::KBase::workspace::ScriptHelpers qw(workspaceURL get_ws_client workspace parseObjectMeta parseWorkspaceMeta);
use Digest::MD5;
$|=1;

my $genomes = {
	471853.12 => [17290,'Beutenbergia cavernae DSM 12333 ','Beutenbergia_cavernae_DSM_12333_'],
	561180.12 => [17290,'Bifidobacterium gallicum DSM 20093','Bifidobacterium_gallicum_DSM_20093'],
	1207543.7 => [17290,'Bifidobacterium bifidum IPLA 20015','Bifidobacterium_bifidum_IPLA_20015'],
	391904.18 => [17290,'Bifidobacterium longum subsp. infantis ATCC 15697','Bifidobacterium_longum_subsp_infantis_ATCC_15697'],
	518634.16 => [17290,'Bifidobacterium breve DSM 20213','Bifidobacterium_breve_DSM_20213'],
	473819.17 => [17290,'Bifidobacterium dentium ATCC 27678','Bifidobacterium_dentium_ATCC_27678'],
	411481.14 => [17290,'Bifidobacterium adolescentis L2-32','Bifidobacterium_adolescentis_L2_32'],
	566552.15 => [17290,'Bifidobacterium catenulatum DSM 16992','Bifidobacterium_catenulatum_DSM_16992'],
	547043.15 => [17290,'Bifidobacterium pseudocatenulatum DSM 20438','Bifidobacterium_pseudocatenulatum_DSM_20438'],
	6666666.81058 => [17195,'Trueperella pyogenes','Trueperella_pyogenes'],
	649743.7 => [17195,'Actinomyces sp. oral taxon 848 str. F0332','Actinomyces_sp_oral_taxon_848_str_F0332'],
	883067.5 => [17195,'Actinobaculum schaalii ','Actinobaculum_schaalii_'],
	883066.5 => [17195,'Actinobaculum massiliae ACS-171-V-Col2','Actinobaculum_massiliae_ACS_171_V_Col2'],
	1321817.5 => [17196,'Actinomyces graevenitzii F0530','Actinomyces_graevenitzii_F0530'],
	521393.6 => [17196,'Actinomyces timonensis DSM 23838','Actinomyces_timonensis_DSM_23838'],
	525246.7 => [17196,'Actinomyces urogenitalis DSM 15434','Actinomyces_urogenitalis_DSM_15434'],
	1120941.6 => [17196,'Actinomyces dentalis DSM 19115','Actinomyces_dentalis_DSM_19115'],
	888056.7 => [17196,'Actinomyces sp. oral taxon 448 str. F0400','Actinomyces_sp_oral_taxon_448_str_F0400'],
	1120943.5 => [17196,'Actinomyces gerencseriae DSM 6844','Actinomyces_gerencseriae_DSM_6844'],
	1167628.6 => [17196,'Actinomyces massiliensis 4401292','Actinomyces_massiliensis_4401292'],
	1278298.5 => [17196,'Actinomyces slackii','Actinomyces_slackii'],
	6666666.81431 => [17196,'Actinomyces oris MG-1','Actinomyces_oris_MG_1'],
	706439.7 => [17196,'Actinomyces sp. oral taxon 171 str. F0337','Actinomyces_sp_oral_taxon_171_str_F0337'],
	653386.7 => [17196,'Actinomyces sp. oral taxon 849 str. F0330','Actinomyces_sp_oral_taxon_849_str_F0330'],
	762963.7 => [17196,'Actinomyces sp. oral taxon 170 str. F0386','Actinomyces_sp_oral_taxon_170_str_F0386'],
	585198.6 => [17197,'Mobiluncus curtisii subsp. curtisii ATCC 35241','Mobiluncus_curtisii_subsp_curtisii_ATCC_35241'],
	871571.6 => [17197,'Mobiluncus mulieris ATCC 35239','Mobiluncus_mulieris_ATCC_35239'],
	1120945.5 => [17197,'Actinomyces neuii BVS029A5','Actinomyces_neuii_BVS029A5'],
	1123488.6 => [17197,'Varibaculum cambriense DSM 15806','Varibaculum_cambriense_DSM_15806'],
	525245.7 => [17198,'Actinomyces coleocanis','Actinomyces_coleocanis'],
	883069.5 => [17198,'Actinomyces europaeus','Actinomyces_europaeus'],
	1120946.5 => [17198,'Actinomyces suimastitidis DSM 15538','Actinomyces_suimastitidis_DSM_15538'],
	1118058.5 => [17198,'Actinomyces sp. ph3','Actinomyces_sp_ph3'],
	1127690.6 => [17198,'Actinomyces sp. oral taxon 181 str. F0379','Actinomyces_sp_oral_taxon_181_str_F0379'],
	1120947.6 => [17198,'Actinomyces vaccimaxillae','Actinomyces_vaccimaxillae'],
	883077.6 => [17198,'Actinomyces turicensis ACS-279-V-Col4','Actinomyces_turicensis_ACS_279_V_Col4'],
	888050.7 => [17198,'Actinomyces cardiffensis F0333','Actinomyces_cardiffensis_F0333'],
	1125717.6 => [17198,'Actinomyces georgiae F0490','Actinomyces_georgiae_F0490'],
	888051.7 => [17198,'Actinomyces sp. oral taxon 178 str. F0338','Actinomyces_sp_oral_taxon_178_str_F0338'],
	888052.6 => [17198,'Actinomyces sp. oral taxon 180 str. F0310','Actinomyces_sp_oral_taxon_180_str_F0310'],
	1105029.6 => [17198,'Actinomyces sp. ICM39','Actinomyces_sp_ICM39'],
	411466.12 => [17198,'Actinomyces odontolyticus ATCC 17982','Actinomyces_odontolyticus_ATCC_17982']
};

my $wsClient = get_ws_client();
my $mssvr = Bio::ModelSEED::MSSeedSupportServer::MSSeedSupportClient->new("http://bio-data-1.mcs.anl.gov/services/ms_fba");
foreach my $genome (keys(%{$genomes})) {
	print "Now loading ".$genome."\n";
	my $workspace = $genomes->{$genome}->[0];
	my $data = $mssvr->getRastGenomeData({
		genome => $genome,
		username => "chenry",
		password => "hello824",
		getSequences => 1,
		getDNASequence => 1
	});
	if (!defined($data->{owner})) {
		die "RAST genome ".$genome." not found!",'get_genomeobject';
	}
	my $genomeObj = {
		id => $genomes->{$genome}->[2],
		scientific_name => $genomes->{$genome}->[1],
		domain => $data->{taxonomy},
		genetic_code => 11,
		dna_size => $data->{size},
		num_contigs => 0,
		contig_lengths => [],
		contig_ids => [],
		source => "RAST",
		source_id => $genome,
		taxonomy => $data->{taxonomy},
		gc_content => 0.5,
		complete => 1,
		publications => [],
		features => [],
    };
    my $contigset = {
		name => $genomeObj->{scientific_name},
		source_id => $genomeObj->{source_id},
		source => $genomeObj->{source},
		type => "Organism",
		contigs => []
    };
    my $contighash = {};
	for (my $i=0; $i < @{$data->{features}}; $i++) {
		my $ftr = $data->{features}->[$i];
		my $feature = {
  			id => $ftr->{ID}->[0],
			type => "peg",
			publications => [],
			subsystems => [],
			protein_families => [],
			aliases => [],
			annotations => [],
			subsystem_data => [],
			regulon_data => [],
			atomic_regulons => [],
			coexpressed_fids => [],
			co_occurring_fids => [],
			protein_translation_length => 0,
			protein_translation => "",
			dna_sequence_length => 0,
			md5 => ""
  		};
  		if ($ftr->{ID}->[0] =~ m/\.([^\.]+)\.\d+$/) {
  			$feature->{type} = $1;
  		}
  		if (defined($ftr->{SEQUENCE})) {
			$feature->{protein_translation} = $ftr->{SEQUENCE}->[0];
			$feature->{protein_translation_length} = length($feature->{protein_translation});
  			$feature->{dna_sequence_length} = 3*$feature->{protein_translation_length};
  			$feature->{md5} = Digest::MD5::md5_hex($feature->{protein_translation});
		}
		if (defined($ftr->{ROLES})) {
			$feature->{function} = join(" / ",@{$ftr->{ROLES}});
		}
  		if (defined($ftr->{LOCATION}->[0]) && $ftr->{LOCATION}->[0] =~ m/^(.+)_(\d+)([\+\-_])(\d+)$/) {
			my $contigData = $1;
			if (!defined($contighash->{$contigData})) {
				$contighash->{$contigData} = $2;
			} elsif ($2 > $contighash->{$contigData}) {
				$contighash->{$contigData} = $2;
			}
			if ($3 eq "-" || $3 eq "+") {
				$feature->{location} = [[$contigData,$2,$3,$4]];
			} elsif ($2 > $4) {
				$feature->{location} = [[$contigData,$2,"-",($2-$4)]];
			} else {
				$feature->{location} = [[$contigData,$2,"+",($4-$2)]];
			}
			$feature->{location}->[0]->[1] = $feature->{location}->[0]->[1]+0;
			$feature->{location}->[0]->[3] = $feature->{location}->[0]->[3]+0;
		}
  		push(@{$genomeObj->{features}},$feature);
	}
	my $ContigObj;
	if (defined($data->{DNAsequence}->[0])) {
    	my $gccount = 0;
    	my $size = 0;
    	for (my $i=0; $i < @{$data->{DNAsequence}}; $i++) {
    		my $closest;
    		foreach my $key (keys(%{$contighash})) {
    			my $dist = abs(length($data->{DNAsequence}->[$i]) - $contighash->{$key});
    			my $closestdist = abs(length($data->{DNAsequence}->[$i]) - $contighash->{$closest});
    			if (!defined($closest) || $dist < $closestdist) {
    				$closest = $key;
    			}
    		}
    		push(@{$contigset->{contigs}},{
    			id => $closest,
				"length" => length($data->{DNAsequence}->[$i]),
				md5 => Digest::MD5::md5_hex($data->{DNAsequence}->[$i]),
				sequence => $data->{DNAsequence}->[$i],
				name => $closest
    		});
    		push(@{$genomeObj->{contig_lengths}},length($data->{DNAsequence}->[$i]));
    		$size += length($data->{DNAsequence}->[$i]);
    		push(@{$genomeObj->{contig_ids}},$closest);
			for ( my $j = 0 ; $j < length($data->{DNAsequence}->[$i]) ; $j++ ) {
				if ( substr( $data->{DNAsequence}->[$i], $j, 1 ) =~ m/[gcGC]/ ) {
					$gccount++;
				}
			}
    	}
    	if ($size > 0) {
			$genomeObj->{gc_content} = $$gccount/$size;
		}
		my $sortedcontigs = [sort { $a->{sequence} cmp $b->{sequence} } @{$contigset->{contigs}}];
		my $str = "";
		for (my $i=0; $i < @{$sortedcontigs}; $i++) {
			if (length($str) > 0) {
				$str .= ";";
			}
			$str .= $sortedcontigs->[$i]->{sequence};	
		}
		$genomeObj->{md5} = Digest::MD5::md5_hex($str);
		$contigset->{md5} = $genomeObj->{md5};
		$contigset->{id} = $genomes->{$genome}->[2].".contigs";
		my $metadata = $wsClient->save_objects({
			'id' => $workspace,
			'objects' => [{
			    type => 'KBaseGenomes.ContigSet',
			    data => $contigset,
			    name => $genomes->{$genome}->[2].".contigs",
			    'meta' => {},
			    'provenance' => []
			}]
		});
		$genomeObj->{contigset_ref} = $metadata->[0]->[6]."/".$metadata->[0]->[0]."/".$metadata->[0]->[4];
	}
	my $metadata = $wsClient->save_objects({
		'id' => $workspace,
		'objects' => [{
		    type => 'KBaseGenomes.Genome',
		    data => $genomeObj,
		    name => $genomes->{$genome}->[2],
		    'meta' => {},
		    'provenance' => []
		}]
	});
}

1;