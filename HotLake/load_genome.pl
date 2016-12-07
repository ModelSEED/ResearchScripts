#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(getUser printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

if (!defined($ARGV[3])) {
	print "Insufficient arguments provided. Usage:\nperl load_genome.pl <fasta file> <gff file> <workspace id> <genome name>\n";
}

my $fasta = $ARGV[0];
my $gff_file = $ARGV[1];
my $workspace = $ARGV[2];
my $genome_name = $ARGV[3];
my $object_id = $ARGV[4];

my @aa_1_letter_order = qw( A C D E F G H I K L M N P Q R S T V W Y );  # Alpha by 1 letter
my @aa_3_letter_order = qw( A R N D C Q E G H I L K M F P S T W Y V );  # PAM matrix order
my @aa_n_codon_order  = qw( L R S A G P T V I C D E F H K N Q Y M W );  

my %genetic_code = (
    # DNA version
    TTT => 'F',  TCT => 'S',  TAT => 'Y',  TGT => 'C',
    TTC => 'F',  TCC => 'S',  TAC => 'Y',  TGC => 'C',
    TTA => 'L',  TCA => 'S',  TAA => '*',  TGA => '*',
    TTG => 'L',  TCG => 'S',  TAG => '*',  TGG => 'W',
    CTT => 'L',  CCT => 'P',  CAT => 'H',  CGT => 'R',
    CTC => 'L',  CCC => 'P',  CAC => 'H',  CGC => 'R',
    CTA => 'L',  CCA => 'P',  CAA => 'Q',  CGA => 'R',
    CTG => 'L',  CCG => 'P',  CAG => 'Q',  CGG => 'R',
    ATT => 'I',  ACT => 'T',  AAT => 'N',  AGT => 'S',
    ATC => 'I',  ACC => 'T',  AAC => 'N',  AGC => 'S',
    ATA => 'I',  ACA => 'T',  AAA => 'K',  AGA => 'R',
    ATG => 'M',  ACG => 'T',  AAG => 'K',  AGG => 'R',
    GTT => 'V',  GCT => 'A',  GAT => 'D',  GGT => 'G',
    GTC => 'V',  GCC => 'A',  GAC => 'D',  GGC => 'G',
    GTA => 'V',  GCA => 'A',  GAA => 'E',  GGA => 'G',
    GTG => 'V',  GCG => 'A',  GAG => 'E',  GGG => 'G',

    #  The following ambiguous encodings are not necessary,  but
    #  speed up the processing of some ambiguous triplets:

    TTY => 'F',  TCY => 'S',  TAY => 'Y',  TGY => 'C',
    TTR => 'L',  TCR => 'S',  TAR => '*',
                 TCN => 'S',
    CTY => 'L',  CCY => 'P',  CAY => 'H',  CGY => 'R',
    CTR => 'L',  CCR => 'P',  CAR => 'Q',  CGR => 'R',
    CTN => 'L',  CCN => 'P',               CGN => 'R',
    ATY => 'I',  ACY => 'T',  AAY => 'N',  AGY => 'S',
                 ACR => 'T',  AAR => 'K',  AGR => 'R',
                 ACN => 'T',
    GTY => 'V',  GCY => 'A',  GAY => 'D',  GGY => 'G',
    GTR => 'V',  GCR => 'A',  GAR => 'E',  GGR => 'G',
    GTN => 'V',  GCN => 'A',               GGN => 'G'
);

#  Add RNA by construction:

foreach ( grep { /T/ } keys %genetic_code )
{
    my $codon = $_;
    $codon =~ s/T/U/g;
    $genetic_code{ $codon } = lc $genetic_code{ $_ }
}

#  Add lower case by construction:

foreach ( keys %genetic_code )
{
    $genetic_code{ lc $_ } = lc $genetic_code{ $_ }
}

#  Construct the genetic code with selenocysteine by difference:

my %genetic_code_with_U = %genetic_code;
$genetic_code_with_U{ TGA } = 'U';
$genetic_code_with_U{ tga } = 'u';
$genetic_code_with_U{ UGA } = 'U';
$genetic_code_with_U{ uga } = 'u';

my %amino_acid_codons_DNA = (
         L  => [ qw( TTA TTG CTA CTG CTT CTC ) ],
         R  => [ qw( AGA AGG CGA CGG CGT CGC ) ],
         S  => [ qw( AGT AGC TCA TCG TCT TCC ) ],
         A  => [ qw( GCA GCG GCT GCC ) ],
         G  => [ qw( GGA GGG GGT GGC ) ],
         P  => [ qw( CCA CCG CCT CCC ) ],
         T  => [ qw( ACA ACG ACT ACC ) ],
         V  => [ qw( GTA GTG GTT GTC ) ],
         I  => [ qw( ATA ATT ATC ) ],
         C  => [ qw( TGT TGC ) ],
         D  => [ qw( GAT GAC ) ],
         E  => [ qw( GAA GAG ) ],
         F  => [ qw( TTT TTC ) ],
         H  => [ qw( CAT CAC ) ],
         K  => [ qw( AAA AAG ) ],
         N  => [ qw( AAT AAC ) ],
         Q  => [ qw( CAA CAG ) ],
         Y  => [ qw( TAT TAC ) ],
         M  => [ qw( ATG ) ],
         U  => [ qw( TGA ) ],
         W  => [ qw( TGG ) ],
         l  => [ qw( tta ttg cta ctg ctt ctc ) ],
         r  => [ qw( aga agg cga cgg cgt cgc ) ],
         s  => [ qw( agt agc tca tcg tct tcc ) ],
         a  => [ qw( gca gcg gct gcc ) ],
         g  => [ qw( gga ggg ggt ggc ) ],
         p  => [ qw( cca ccg cct ccc ) ],
         t  => [ qw( aca acg act acc ) ],
         v  => [ qw( gta gtg gtt gtc ) ],
         i  => [ qw( ata att atc ) ],
         c  => [ qw( tgt tgc ) ],
         d  => [ qw( gat gac ) ],
         e  => [ qw( gaa gag ) ],
         f  => [ qw( ttt ttc ) ],
         h  => [ qw( cat cac ) ],
         k  => [ qw( aaa aag ) ],
         n  => [ qw( aat aac ) ],
         q  => [ qw( caa cag ) ],
         y  => [ qw( tat tac ) ],
         m  => [ qw( atg ) ],
         u  => [ qw( tga ) ],
         w  => [ qw( tgg ) ],
        '*' => [ qw( TAA TAG TGA ) ]
);

my %amino_acid_codons_RNA = (
         L  => [ qw( UUA UUG CUA CUG CUU CUC ) ],
         R  => [ qw( AGA AGG CGA CGG CGU CGC ) ],
         S  => [ qw( AGU AGC UCA UCG UCU UCC ) ],
         A  => [ qw( GCA GCG GCU GCC ) ],
         G  => [ qw( GGA GGG GGU GGC ) ],
         P  => [ qw( CCA CCG CCU CCC ) ],
         T  => [ qw( ACA ACG ACU ACC ) ],
         V  => [ qw( GUA GUG GUU GUC ) ],
         B  => [ qw( GAU GAC AAU AAC ) ],
         Z  => [ qw( GAA GAG CAA CAG ) ],
         I  => [ qw( AUA AUU AUC ) ],
         C  => [ qw( UGU UGC ) ],
         D  => [ qw( GAU GAC ) ],
         E  => [ qw( GAA GAG ) ],
         F  => [ qw( UUU UUC ) ],
         H  => [ qw( CAU CAC ) ],
         K  => [ qw( AAA AAG ) ],
         N  => [ qw( AAU AAC ) ],
         Q  => [ qw( CAA CAG ) ],
         Y  => [ qw( UAU UAC ) ],
         M  => [ qw( AUG ) ],
         U  => [ qw( UGA ) ],
         W  => [ qw( UGG ) ],
         l  => [ qw( uua uug cua cug cuu cuc ) ],
         r  => [ qw( aga agg cga cgg cgu cgc ) ],
         s  => [ qw( agu agc uca ucg ucu ucc ) ],
         a  => [ qw( gca gcg gcu gcc ) ],
         g  => [ qw( gga ggg ggu ggc ) ],
         p  => [ qw( cca ccg ccu ccc ) ],
         t  => [ qw( aca acg acu acc ) ],
         v  => [ qw( gua gug guu guc ) ],
         b  => [ qw( gau gac aau aac ) ],
         z  => [ qw( gaa gag caa cag ) ],
         i  => [ qw( aua auu auc ) ],
         c  => [ qw( ugu ugc ) ],
         d  => [ qw( gau gac ) ],
         e  => [ qw( gaa gag ) ],
         f  => [ qw( uuu uuc ) ],
         h  => [ qw( cau cac ) ],
         k  => [ qw( aaa aag ) ],
         n  => [ qw( aau aac ) ],
         q  => [ qw( caa cag ) ],
         y  => [ qw( uau uac ) ],
         m  => [ qw( aug ) ],
         u  => [ qw( uga ) ],
         w  => [ qw( ugg ) ],
        '*' => [ qw( UAA UAG UGA ) ]
);


my %n_codon_for_aa = map {
    $_ => scalar @{ $amino_acid_codons_DNA{ $_ } }
    } keys %amino_acid_codons_DNA;

my %reverse_genetic_code_DNA = (
         A  => "GCN",  a  => "gcn",
         C  => "TGY",  c  => "tgy",
         D  => "GAY",  d  => "gay",
         E  => "GAR",  e  => "gar",
         F  => "TTY",  f  => "tty",
         G  => "GGN",  g  => "ggn",
         H  => "CAY",  h  => "cay",
         I  => "ATH",  i  => "ath",
         K  => "AAR",  k  => "aar",
         L  => "YTN",  l  => "ytn",
         M  => "ATG",  m  => "atg",
         N  => "AAY",  n  => "aay",
         P  => "CCN",  p  => "ccn",
         Q  => "CAR",  q  => "car",
         R  => "MGN",  r  => "mgn",
         S  => "WSN",  s  => "wsn",
         T  => "ACN",  t  => "acn",
         U  => "TGA",  u  => "tga",
         V  => "GTN",  v  => "gtn",
         W  => "TGG",  w  => "tgg",
         X  => "NNN",  x  => "nnn",
         Y  => "TAY",  y  => "tay",
        '*' => "TRR"
);

my %reverse_genetic_code_RNA = (
         A  => "GCN",  a  => "gcn",
         C  => "UGY",  c  => "ugy",
         D  => "GAY",  d  => "gay",
         E  => "GAR",  e  => "gar",
         F  => "UUY",  f  => "uuy",
         G  => "GGN",  g  => "ggn",
         H  => "CAY",  h  => "cay",
         I  => "AUH",  i  => "auh",
         K  => "AAR",  k  => "aar",
         L  => "YUN",  l  => "yun",
         M  => "AUG",  m  => "aug",
         N  => "AAY",  n  => "aay",
         P  => "CCN",  p  => "ccn",
         Q  => "CAR",  q  => "car",
         R  => "MGN",  r  => "mgn",
         S  => "WSN",  s  => "wsn",
         T  => "ACN",  t  => "acn",
         U  => "UGA",  u  => "uga",
         V  => "GUN",  v  => "gun",
         W  => "UGG",  w  => "ugg",
         X  => "NNN",  x  => "nnn",
         Y  => "UAY",  y  => "uay",
        '*' => "URR"
);

my %DNA_letter_can_be = (
    A => ["A"],                 a => ["a"],
    B => ["C", "G", "T"],       b => ["c", "g", "t"],
    C => ["C"],                 c => ["c"],
    D => ["A", "G", "T"],       d => ["a", "g", "t"],
    G => ["G"],                 g => ["g"],
    H => ["A", "C", "T"],       h => ["a", "c", "t"],
    K => ["G", "T"],            k => ["g", "t"],
    M => ["A", "C"],            m => ["a", "c"],
    N => ["A", "C", "G", "T"],  n => ["a", "c", "g", "t"],
    R => ["A", "G"],            r => ["a", "g"],
    S => ["C", "G"],            s => ["c", "g"],
    T => ["T"],                 t => ["t"],
    U => ["T"],                 u => ["t"],
    V => ["A", "C", "G"],       v => ["a", "c", "g"],
    W => ["A", "T"],            w => ["a", "t"],
    Y => ["C", "T"],            y => ["c", "t"]
);

my %RNA_letter_can_be = (
    A => ["A"],                 a => ["a"],
    B => ["C", "G", "U"],       b => ["c", "g", "u"],
    C => ["C"],                 c => ["c"],
    D => ["A", "G", "U"],       d => ["a", "g", "u"],
    G => ["G"],                 g => ["g"],
    H => ["A", "C", "U"],       h => ["a", "c", "u"],
    K => ["G", "U"],            k => ["g", "u"],
    M => ["A", "C"],            m => ["a", "c"],
    N => ["A", "C", "G", "U"],  n => ["a", "c", "g", "u"],
    R => ["A", "G"],            r => ["a", "g"],
    S => ["C", "G"],            s => ["c", "g"],
    T => ["U"],                 t => ["u"],
    U => ["U"],                 u => ["u"],
    V => ["A", "C", "G"],       v => ["a", "c", "g"],
    W => ["A", "U"],            w => ["a", "u"],
    Y => ["C", "U"],            y => ["c", "u"]
);

my %one_letter_to_three_letter_aa = (
         A  => "Ala", a  => "Ala",
         B  => "Asx", b  => "Asx",
         C  => "Cys", c  => "Cys",
         D  => "Asp", d  => "Asp",
         E  => "Glu", e  => "Glu",
         F  => "Phe", f  => "Phe",
         G  => "Gly", g  => "Gly",
         H  => "His", h  => "His",
         I  => "Ile", i  => "Ile",
         K  => "Lys", k  => "Lys",
         L  => "Leu", l  => "Leu",
         M  => "Met", m  => "Met",
         N  => "Asn", n  => "Asn",
         P  => "Pro", p  => "Pro",
         Q  => "Gln", q  => "Gln",
         R  => "Arg", r  => "Arg",
         S  => "Ser", s  => "Ser",
         T  => "Thr", t  => "Thr",
         U  => "Sec", u  => "Sec",
         V  => "Val", v  => "Val",
         W  => "Trp", w  => "Trp",
         X  => "Xxx", x  => "Xxx",
         Y  => "Tyr", y  => "Tyr",
         Z  => "Glx", z  => "Glx",
        '*' => "***"
        );

my %three_letter_to_one_letter_aa = (
     ALA  => "A",   Ala  => "A",   ala  => "a",
     ARG  => "R",   Arg  => "R",   arg  => "r",
     ASN  => "N",   Asn  => "N",   asn  => "n",
     ASP  => "D",   Asp  => "D",   asp  => "d",
     ASX  => "B",   Asx  => "B",   asx  => "b",
     CYS  => "C",   Cys  => "C",   cys  => "c",
     GLN  => "Q",   Gln  => "Q",   gln  => "q",
     GLU  => "E",   Glu  => "E",   glu  => "e",
     GLX  => "Z",   Glx  => "Z",   glx  => "z",
     GLY  => "G",   Gly  => "G",   gly  => "g",
     HIS  => "H",   His  => "H",   his  => "h",
     ILE  => "I",   Ile  => "I",   ile  => "i",
     LEU  => "L",   Leu  => "L",   leu  => "l",
     LYS  => "K",   Lys  => "K",   lys  => "k",
     MET  => "M",   Met  => "M",   met  => "m",
     PHE  => "F",   Phe  => "F",   phe  => "f",
     PRO  => "P",   Pro  => "P",   pro  => "p",
     SEC  => "U",   Sec  => "U",   sec  => "u",
     SER  => "S",   Ser  => "S",   ser  => "s",
     THR  => "T",   Thr  => "T",   thr  => "t",
     TRP  => "W",   Trp  => "W",   trp  => "w",
     TYR  => "Y",   Tyr  => "Y",   tyr  => "y",
     VAL  => "V",   Val  => "V",   val  => "v",
     XAA  => "X",   Xaa  => "X",   xaa  => "x",
     XXX  => "X",   Xxx  => "X",   xxx  => "x",
    '***' => "*"
);

#Loading FASTA object
my $object = {
	id => $object_id.".contigs",
	name => $genome_name,
	source_id => $object_id.".contigs",
	source => getUser(),
	type => "Organism",
	contigs => [],
	md5 => undef
};
open FILE, "<",$fasta  or do {
	die "$0: open ".$fasta.": $!";
};
my $contig;
my $sequence = "";
my $description;
my $length = 0;
my $lengths = [];
my $contigids = [];
my $contigseq = {};
my $gc = 0;
my $numcontigs = 0;
while (my $line = <FILE>) {
	chomp($line);
	if ($line =~ m/\>(.+)/) {
		if (defined($contig)) {
			$numcontigs++;
			push(@{$lengths},length($sequence));
			push(@{$contigids},$contig);
			$length += length($sequence);
			$contigseq->{$contig} = $sequence;
			push(@{$object->{contigs}},{
				id => $contig,
				"length" => length($sequence),
				md5 => Digest::MD5::md5_hex($sequence),
				sequence => $sequence,
				genetic_code => 11,
				name => $contig,
				complete => 1,
				description => $description
			});
			$sequence =~ s/[atAT]//g;
			$gc += length($sequence);
		}
		$sequence = "";
		$contig = $1;
		$description = $1;
		if($contig =~ /(.*)?\s(.+)/ ) {
			$contig = $1;
			$description = $2;
		}
	} else {
		$sequence .= $line;
	}
}
if (defined($contig)) {
	$numcontigs++;
	push(@{$lengths},length($sequence));
	push(@{$contigids},$contig);
	$length += length($sequence);
	$contigseq->{$contig} = $sequence;
	push(@{$object->{contigs}},{
		id => $contig,
		"length" => length($sequence),
		md5 => Digest::MD5::md5_hex($sequence),
		sequence => $sequence,
		genetic_code => 11,
		name => $contig,
		complete => 1,
		description => $description
	});
	$sequence =~ s/[atAT]//g;
	$gc += length($sequence);
}
$gc = $gc/$length;
close(FILE);
$object->{contigs} = [sort { $a->{sequence} cmp $b->{sequence} } @{$object->{contigs}}];
my $str = "";
for (my $i=0; $i < @{$object->{contigs}}; $i++) {
	if (length($str) > 0) {
		$str .= ";";
	}
	$str .= $object->{contigs}->[$i]->{sequence};
}
$object->{md5} = Digest::MD5::md5_hex($str);
save_workspace_object($workspace."/".$object->{id},$object,"KBaseGenomes.ContigSet");
#Loading Genome object
my $md5 = $object->{md5};
$object = {
	id => $object_id,
	scientific_name => $genome_name,
	domain => "Bacteria",
	genetic_code => 11,
	dna_size => $length,
	num_contigs => $numcontigs,
	contig_lengths => $lengths,
	contig_ids => $contigids,
	source => "HotLake SFA",
	source_id => $object_id,
	md5 => $object->{md5},
	taxonomy => "Unknown",
	gc_content => $gc,
	complete => 1,
	features => [],
	contigset_ref => $workspace."/".$object->{id}
};
open(FILEF, "<", $gff_file);
while (my $line = <FILEF>) {
	chomp($line);
	if (defined($line)) {
		my $array = [split(/\t/,$line)];
		my $ftr = {
			location => [[$array->[0],$array->[3],$array->[6],($array->[4]-$array->[3])]],
			type => "CDS",
			function => "Unknown",
			aliases => [],
			annotations => []
		};
		if ($line =~ m/orig_coding_type\s\"([^\"]+)\";/) {
			$ftr->{type} = $1;
		}
		if ($line =~ m/gene_id\s\"([^\"]+)\";/) {
			$ftr->{id} = $1;
		}
		if ($line =~ m/product\s\"([^\"]+)\";/) {
			$ftr->{function} = $1;
		}
		if ($line =~ m/EC_number\s\"([^\"]+)\";/) {
			$ftr->{function} .= " [".$1."]";
		}
		if ($line =~ m/transcript_id\s\"([^\"]+)\";/) {
			push(@{$ftr->{aliases}},$1);
		}
		if ($line =~ m/locus_tag\s\"([^\"]+)\";/) {
			$ftr->{id} = $1;
			#push(@{$ftr->{aliases}},$1);
		}
		if ($line =~ m/RAST\s\"([^\"]+)\";/) {
			push(@{$ftr->{aliases}},$1);
		}
		my $seq;
		if ($array->[6] eq "-") {
			$seq = scalar reverse substr($contigseq->{$array->[0]},$array->[3],($array->[4]-$array->[3]));
			$seq =~ s/A/M/g;
			$seq =~ s/a/m/g;
			$seq =~ s/T/A/g;
			$seq =~ s/t/a/g;
			$seq =~ s/M/T/g;
			$seq =~ s/m/t/g;
			$seq =~ s/G/M/g;
			$seq =~ s/g/m/g;
			$seq =~ s/C/G/g;
			$seq =~ s/c/g/g;
			$seq =~ s/M/C/g;
			$seq =~ s/m/c/g;
			$ftr->{location}->[0]->[1] = $array->[4];
		} else {
			$seq = substr($contigseq->{$array->[0]},$array->[3]-1,($array->[4]-$array->[3]));
		}
		$ftr->{location}->[0]->[1] = $ftr->{location}->[0]->[1]+0;
		$ftr->{location}->[0]->[3] = $ftr->{location}->[0]->[3]+0;
		$ftr->{dna_sequence} = $seq;
		$ftr->{dna_sequence_length} = length($ftr->{dna_sequence});
		my $proteinSeq = &translate_seq($ftr->{dna_sequence});
		$proteinSeq =~ s/[xX]+$//;
		$ftr->{protein_translation} = $proteinSeq;
		$ftr->{protein_translation_length} = length($ftr->{protein_translation});
		$ftr->{md5} = Digest::MD5::md5_hex($ftr->{protein_translation});
		push(@{$object->{features}},$ftr);
	}
}
close(FILEF);
save_workspace_object($workspace."/".$object->{id},$object,"KBaseGenomes.Genome");

sub translate_seq {
    my($seq) = @_;
    $seq =~ tr/-//d;        #  remove gaps
    my @codons = $seq =~ m/(...?)/g;  #  Will try to translate last 2 nt
    #  A second argument that is true forces first amino acid to be Met

    my @met;
    if ( ( shift @_ ) && ( my $codon1 = shift @codons ) )
    {
        push @met, ( $codon1 =~ /[a-z]/ ? 'm' : 'M' );
    }

    join( '', @met, map { translate_codon( $_ ) } @codons )
}

sub translate_codon {
    my $codon = shift;
    $codon =~ tr/Uu/Tt/;     #  Make it DNA

    #  Try a simple lookup:

    my $aa;
    if ( $aa = $genetic_code{ $codon } ) { return $aa }

    #  Attempt to recover from mixed-case codons:

    $codon = ( $codon =~ /[a-z]/ ) ? lc $codon : uc $codon;
    if ( $aa = $genetic_code{ $codon } ) { return $aa }

    #  The code defined above catches simple N, R and Y ambiguities in the
    #  third position.  Other codons (e.g., GG[KMSWBDHV], or even GG) might
    #  be unambiguously translated by converting the last position to N and
    #  seeing if this is in the code table:

    my $N = ( $codon =~ /[a-z]/ ) ? 'n' : 'N';
    if ( $aa = $genetic_code{ substr($codon,0,2) . $N } ) { return $aa }

    #  Test that codon is valid for an unambiguous aa:

    my $X = ( $codon =~ /[a-z]/ ) ? 'x' : 'X';
    if ( $codon !~ m/^[ACGTMY][ACGT][ACGTKMRSWYBDHVN]$/i
      && $codon !~ m/^YT[AGR]$/i     #  Leu YTR
      && $codon !~ m/^MG[AGR]$/i     #  Arg MGR
       )
    {
        return $X;
    }

    #  Expand all ambiguous nucleotides to see if they all yield same aa.

    my @n1 = @{ $DNA_letter_can_be{ substr( $codon, 0, 1 ) } };
    my $n2 =                        substr( $codon, 1, 1 );
    my @n3 = @{ $DNA_letter_can_be{ substr( $codon, 2, 1 ) } };
    my @triples = map { my $n12 = $_ . $n2; map { $n12 . $_ } @n3 } @n1;

    my $triple = shift @triples;
    $aa = $genetic_code{ $triple };
    $aa or return $X;

    foreach $triple ( @triples ) { return $X if $aa ne $genetic_code{$triple} }
    return $aa;
}