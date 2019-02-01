#!/usr/bin/perl -w

use strict;
#use JSON::XS;
use JSON;
use POSIX;
use Digest::MD5;
local $| = 1;

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


my $path = $ARGV[0];
my $size = $ARGV[1];
my $genomes = &FROMJSON(join("\n",@{&LOADFILE($path."/AllPatricGenomes.txt")}));
my $totalgenomes = keys(%{$genomes});

my $count = 0;
my $tooshort = 0;
my $overlapping_sketches = {
	"000" => 0,
	"100" => 0,
	"010" => 0,
	"001" => 0,
	"110" => 0,
	"011" => 0,
	"101" => 0,
	"111" => 0
};
my $hash = {};
my $funchash = {};
my $genomelist = [keys(%{$genomes})];
print "Sketch start:".time()."\n";
for (my $i=0; $i < @{$genomelist}; $i++) { 
	print "Sketching:".$i." of ".$totalgenomes."\t".time()."\n";
	my $lines = &LOADFILE("/vol/patric3/downloads/genomes/".$genomelist->[$i]."/".$genomelist->[$i].".PATRIC.faa");
	my $id;
	my $func;
	my $seq;
	for (my $j=0; $j < @{$lines}; $j++) {
		if ($lines->[$j] =~ m/^>([^\s^\t]+)/) {
			my $newid = $1;
			if (defined($id)) {
				if (50+$size < length($seq)) {
					my $overlap = "";
					my $sketch = substr($seq,49,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap = "1";
					} else {
						$overlap = "0";
					}
					push(@{$hash->{$sketch}},[$id,0]);
					$sketch = substr($seq,length($seq)-49-$size,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap .= "1";
					} else {
						$overlap .= "0";
					}
					push(@{$hash->{$sketch}},[$id,1]);
					$sketch = substr($seq,length($seq)/2-$size/2,$size);
					if ($size > 30) {
						$sketch = Digest::MD5::md5_hex($sketch);
					}
					if (defined($hash->{$sketch})) {
						$overlap .= "1";
					} else {
						$overlap .= "0";
					}
					$overlapping_sketches->{$overlap}++;
					push(@{$hash->{$sketch}},[$id,2]);
				} else {
					$tooshort++;
				}
				$count++;
			}
			if ($lines->[$j] =~ m/^>([^\s^\t]+)\s\s\s(.+)\s\s\s/) {
				$funchash->{$newid} = $2;
			}
			$id = $newid;
			$seq = "";
		} else {
			$seq .= $lines->[$j];
		}
	}
#	if ($count > 1000) {
#		last;
#	}
}
print "Count:".$count."\n";
print "Too short:".$tooshort."\n";
foreach my $key (keys(%{$overlapping_sketches})) {
	print $key."\t".$overlapping_sketches->{$key}."\n";
}

my $filelist = [qw(
GO-all
GO-ERR2162200
GO-ERR2162201
GO-ERR2162202
GO-ERR2162203
GO-ERR2162204
GO-ERR2162205
GO-ERR2162206
GO-ERR2162207
GO-ERR2162208
GO-ERR2162209
GO-ERR2162210
GO-ERR2162211
GO-ERR2162212
GO-ERR2162213
GO-ERR2162214
GO-ERR2162215
GO-ERR2162216
GO-ERR2162217
GO-ERR2162218
GO-ERR2162219
GO-ERR2162220
GO-ERR2162221
GO-ERR2162222
GO-ERR2162223
GO-ERR2162224
)];
my $output = {};
for (my $i=0; $i < @{$filelist}; $i++) {
	print "Gene scan ".$filelist->[$i]."\n";
	my $features = &FROMJSON(join("\n",@{&LOADFILE($path."/all_genes/".$filelist->[$i].".json")}));
	foreach my $contigid (keys(%{$features})) {
		$output->{$filelist->[$i]}->{$contigid} = {
			genes => {},
			raw => {},
			genehits => 0,
			rawhits => 0
		};
		for (my $j=0; $j < @{$features->{$contigid}}; $j++) {
			$output->{$filelist->[$i]}->{$contigid}->{genes}->{$features->{$contigid}->[$j]->{id}} = {
				s => {},
				f => {},
				count => 0
			};
			my $seq = $features->{$contigid}->[$j]->{protein_translation};
			my $high = $features->{$contigid}->[$j]->{location}->[0]->[1] + $features->{$contigid}->[$j]->{location}->[0]->[3];
			my $low = $features->{$contigid}->[$j]->{location}->[0]->[1];
			if ($features->{$contigid}->[$j]->{location}->[0]->[2] eq "-") {
				$low = $features->{$contigid}->[$j]->{location}->[0]->[1] - $features->{$contigid}->[$j]->{location}->[0]->[3];
				$high = $features->{$contigid}->[$j]->{location}->[0]->[1];
			}
			if (!defined($output->{$filelist->[$i]}->{$contigid}->{lowest}) || $output->{$filelist->[$i]}->{$contigid}->{lowest} > $low) {
				$output->{$filelist->[$i]}->{$contigid}->{lowest} = $low;
			}
			if (!defined($output->{$filelist->[$i]}->{$contigid}->{highest}) || $output->{$filelist->[$i]}->{$contigid}->{highest} < $high) {
				$output->{$filelist->[$i]}->{$contigid}->{highest} = $high;
			}
			&ScanProteinForHits($seq,$output->{$filelist->[$i]},"protein",$contigid,$features->{$contigid}->[$j]->{id});
			if ($output->{$filelist->[$i]}->{$contigid}->{genes}->{$features->{$contigid}->[$j]->{id}}->{count} > 0) {
				$output->{$filelist->[$i]}->{$contigid}->{genehits}++;
			}
		}
	}
}
for (my $i=0; $i < @{$filelist}; $i++) {
	print "Contig scan ".$filelist->[$i]."\n";
	my $lines = &LOADFILE($path."/all_genes/".$filelist->[$i].".fasta");
	my $id;
	my $func;
	my $seq = "";
	for (my $j=0; $j < @{$lines}; $j++) {
	#for (my $j=0; $j < 1000; $j++) {
		if ($lines->[$j] =~ m/^>(\w+)[\s\t]/) {
			my $newid = $1;
			if (defined($id)) {
				if (defined($output->{$filelist->[$i]}->{$id})) {
					if ($output->{$filelist->[$i]}->{$id}->{lowest} > 3*$size) {
						$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence(substr($seq,0,$output->{$filelist->[$i]}->{$id}->{lowest}),$output->{$filelist->[$i]},"head",$id);
					}
					if (length($seq) > $output->{$filelist->[$i]}->{$id}->{highest}+3*$size) {
						$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence(substr($seq,$output->{$filelist->[$i]}->{$id}->{highest}),$output->{$filelist->[$i]},"tail",$id);
					}
				} else {
					$output->{$filelist->[$i]}->{$id}->{rawhits} = &CheckContigSequence($seq,$output->{$filelist->[$i]},"entire",$id);
				}
			}
			$id = $newid;
			$seq = "";
		} else {
			$seq .= $lines->[$j];
		}
	}	
}

foreach my $id (keys(%{$output})) {
	my $genespeclines = ["Contig id\tGene id\tSpecies\tCount\tFraction"];
	my $genefunclines = ["Contig id\tGene id\tFunction\tCount\tFraction"];
	my $cmlines = ["Contig id\tGenes\tGeneHits\tRawHits\tHigh\tLow"];
	my $rawspeclines = ["Contig id\tType\tSpecies\tCount\tFraction"];
	my $rawfunclines = ["Contig id\tType\tFunction\tCount\tFraction"];
	foreach my $contigid (keys(%{$output->{$id}})) {
		my $genecount = keys(%{$output->{$id}->{$contigid}->{genes}});
		my $high = "-";
		my $low = "-";
		if (defined($output->{$id}->{$contigid}->{lowest})) {
			$low = $output->{$id}->{$contigid}->{lowest};
		}
		if (defined($output->{$id}->{$contigid}->{highest})) {
			$high = $output->{$id}->{$contigid}->{highest};
		}
		push(@{$cmlines},$contigid."\t".$genecount."\t".$output->{$id}->{$contigid}->{genehits}."\t".$output->{$id}->{$contigid}->{rawhits}."\t".$high."\t".$low);
		foreach my $gene (keys(%{$output->{$id}->{$contigid}->{genes}})) {
			foreach my $species (keys(%{$output->{$id}->{$contigid}->{genes}->{$gene}->{s}})) {
				my $count = $output->{$id}->{$contigid}->{genes}->{$gene}->{s}->{$species};
				my $fraction = $count/$output->{$id}->{$contigid}->{genes}->{$gene}->{count};
				push(@{$genespeclines},$contigid."\t".$gene."\t".$species."\t".$count."\t".$fraction);
			}
			foreach my $function (keys(%{$output->{$id}->{$contigid}->{genes}->{$gene}->{f}})) {
				my $count = $output->{$id}->{$contigid}->{genes}->{$gene}->{f}->{$function};
				my $fraction = $count/$output->{$id}->{$contigid}->{genes}->{$gene}->{count};
				push(@{$genefunclines},$contigid."\t".$gene."\t".$function."\t".$count."\t".$fraction);
			}
		}
		foreach my $gene (keys(%{$output->{$id}->{$contigid}->{raw}})) {
			foreach my $species (keys(%{$output->{$id}->{$contigid}->{raw}->{$gene}->{s}})) {
				my $count = $output->{$id}->{$contigid}->{raw}->{$gene}->{s}->{$species};
				my $fraction = $count/$output->{$id}->{$contigid}->{raw}->{$gene}->{count};
				push(@{$rawspeclines},$contigid."\t".$gene."\t".$species."\t".$count."\t".$fraction);
			}
			foreach my $function (keys(%{$output->{$id}->{$contigid}->{raw}->{$gene}->{f}})) {
				my $count = $output->{$id}->{$contigid}->{raw}->{$gene}->{f}->{$function};
				my $fraction = $count/$output->{$id}->{$contigid}->{raw}->{$gene}->{count};
				push(@{$rawfunclines},$contigid."\t".$gene."\t".$function."\t".$count."\t".$fraction);
			}
		}
	}
	&PRINTFILE($path."/".$size."_".$id.".contig_meta",$cmlines);
	&PRINTFILE($path."/".$size."_".$id.".gene_spec",$genespeclines);
	&PRINTFILE($path."/".$size."_".$id.".gene_func",$genefunclines);
	&PRINTFILE($path."/".$size."_".$id.".contig_spec",$rawspeclines);
	&PRINTFILE($path."/".$size."_".$id.".contig_func",$rawfunclines);
}

sub CheckContigSequence {
    my ($seq,$outhash,$type,$contigid) = @_;
	my $count = 0;
	my $protseq = &translate_sequence($seq,1);
	print "1:".$protseq."\n";
	&ScanProteinForHits($protseq,$outhash,$type."_F1",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F1"}) && $outhash->{$contigid}->{raw}->{$type."_F1"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F1"}->{count};
	}
	$protseq = &translate_sequence(substr($seq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F2",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F2"}) && $outhash->{$contigid}->{raw}->{$type."_F2"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F2"}->{count};
	}
	$protseq = &translate_sequence(substr($seq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_F3",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_F3"}) && $outhash->{$contigid}->{raw}->{$type."_F3"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_F3"}->{count};
	}
	my $revseq = &reverse_sequence($seq);
	$protseq = &translate_sequence($revseq,1);
	&ScanProteinForHits($protseq,$outhash,$type."_R1",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_R1"}) && $outhash->{$contigid}->{raw}->{$type."_R1"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R1"}->{count};
	}
	$protseq = &translate_sequence(substr($revseq,1),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R2",$contigid);
	if (defined($outhash->{$contigid}->{raw}->{$type."_R2"}) && $outhash->{$contigid}->{raw}->{$type."_R2"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R2"}->{count};
	}
	$protseq = &translate_sequence(substr($revseq,2),1);
	&ScanProteinForHits($protseq,$outhash,$type."_R3",$contigid);	
	if (defined($outhash->{$contigid}->{raw}->{$type."_R3"}) && $outhash->{$contigid}->{raw}->{$type."_R3"}->{count} > 0) {
		$count += $outhash->{$contigid}->{raw}->{$type."_R3"}->{count};
	}
	return $count;
}

sub ScanProteinForHits {
	my ($seq,$outhash,$type,$contigid,$id) = @_;
	if ($size < length($seq)) {
		my $genehash = {};
		for (my $k=0; $k < (length($seq)-$size-1); $k++) {
			my $query = substr($seq,$k,$size);
			if ($size > 30) {
				$query = Digest::MD5::md5_hex($query);
			}
			my $pointer = $hash->{$query};
			if (defined($pointer)) {
				for (my $m=0; $m < @{$pointer}; $m++) {
					my $new = 0;
					if (!defined($genehash->{$pointer->[$m]->[0]})) {
						$genehash->{$pointer->[$m]->[0]} = 0;
						$new = 1;
					}
					$genehash->{$pointer->[$m]->[0]}++;
					if ($new == 1) {
						my $function = "unknown";
						if (defined($funchash->{$pointer->[$m]->[0]})) {
							$function = $funchash->{$pointer->[$m]->[0]};
						}
						my $species = "unknown";
						if ($pointer->[$m]->[0] =~ m/\|(\d+\.\d+)\./) {
							$species = $genomes->{$1}->{n};
							my $array = [split(/\s/,$species)];
							if ($array->[1] eq "sp.") {
								$species = $array->[0]." ".$array->[2];
							} else {
								$species = $array->[0]." ".$array->[1];
							}
						}
						my $hash;
						if ($type eq "protein") {
							if (!defined($outhash->{$contigid}->{genes}->{$id})) {
								$outhash->{$contigid}->{genes}->{$id} = {s => {},f => {},count => 0};
							}
							$hash = $outhash->{$contigid}->{genes}->{$id};
						} else {
							if (!defined($outhash->{$contigid}->{raw}->{$type})) {
								$outhash->{$contigid}->{raw}->{$type} = {s => {},f => {},count => 0};
							}
							$hash = $outhash->{$contigid}->{raw}->{$type};
						}
						if (!defined($hash->{s}->{$species})) {
							$hash->{s}->{$species} = 0;
						}
						$hash->{s}->{$species}++;
						if (!defined($hash->{f}->{$function})) {
							$hash->{f}->{$function} = 0;
						}
						$hash->{f}->{$function}++;
						$hash->{count}++;
					}	
				}
			}
		}		
	}
}

sub LOADFILE {
    my ($filename) = @_;
    my $DataArrayRef = [];
    open (my $fh, "<", $filename);
    while (my $Line = <$fh>) {
        $Line =~ s/\r//;
        chomp($Line);
        push(@{$DataArrayRef},$Line);
    }
    close($fh);
    return $DataArrayRef;
}

sub PRINTFILE {
    my ($filename,$arrayRef) = @_;
    open ( my $fh, ">", $filename);
    foreach my $Item (@{$arrayRef}) {
    	print $fh $Item."\n";
    }
    close($fh);
}

sub TOJSON {
    my ($ref,$prettyprint) = @_;
    my $JSON = JSON->new->utf8(1);
    if (defined($prettyprint) && $prettyprint == 1) {
		$JSON->pretty(1);
    }
    return $JSON->encode($ref);
}

sub FROMJSON {
    my ($data) = @_;
    return decode_json $data;
}

sub reverse_sequence {
	my($sequence) = @_;
	$sequence = scalar reverse $sequence;
	$sequence =~ s/A/M/g;
	$sequence =~ s/a/m/g;
	$sequence =~ s/T/A/g;
	$sequence =~ s/t/a/g;
	$sequence =~ s/M/T/g;
	$sequence =~ s/m/t/g;
	$sequence =~ s/G/M/g;
	$sequence =~ s/g/m/g;
	$sequence =~ s/C/G/g;
	$sequence =~ s/c/g/g;
	$sequence =~ s/M/C/g;
	$sequence =~ s/m/c/g;
	return $sequence;
}

sub translate_sequence {
	my($sequence,$clearx) = @_;
	
	$sequence =~ tr/-//d;        #  remove gaps
    my @codons = $sequence =~ m/(...?)/g;  #  Will try to translate last 2 nt
    #  A second argument that is true forces first amino acid to be Met

    my @met;
    if ( ( shift @_ ) && ( my $codon1 = shift @codons ) ) {
        push @met, ( $codon1 =~ /[a-z]/ ? 'm' : 'M' );
    }

    my $protseq = join( '', @met, map { translate_codon( $_ ) } @codons );
	if (defined($clearx) && $clearx == 1) {
		$protseq =~ s/[xX]+$//;
	}
	return $protseq;
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