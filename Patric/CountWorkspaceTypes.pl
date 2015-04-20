use JSON::XS;
use Data::Dumper;
use MongoDB::Connection;

my $config = {
	host => "localhost",
	db_name => "WorkspaceBuild",
	auto_connect => 1,
	auto_reconnect => 1
};
my $conn = MongoDB::Connection->new(%$config);
my $db = $conn->get_database("WorkspaceBuild");
my $types = [qw(
string
genome
unspecified
folder
transcriptomics_experiment
proteomics_experiment
genome_group
feature_group
experiment_group
contigs
genome_annotation_result
reads
job_result
diffexp_input_data
diffexp_input_metadata
genbank_file
feature_protein_fasta
feature_dna_fasta
feature_table
gff
embl
diffexp_experiment
diffexp_expression
diffexp_mapping
diffexp_sample
pdf
xls
xlsx
doc
docx
ppt
pptx
html
gif
jpg
png
txt
csv 
tar_gz
)];
print "Counts since conception:\n";
foreach my $type (@{$types}) {
	print $type."\t".$db->get_collection('objects')->count({
		type => $type
	})."\n";
}
print "Counts since march:\n";
foreach my $type (@{$types}) {
	print $type."\t".$db->get_collection('objects')->count({
		creation_date => { '$gt' => "2015-03" },
		type => $type
	})."\n";
}