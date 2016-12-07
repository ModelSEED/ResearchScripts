use Data::MATFile qw/read_matfile/;
use Data::Dumper;

my $filename = $ARGV[0];

my $matfile = read_matfile ($filename);

print Data::Dumper->Dump($matfile);