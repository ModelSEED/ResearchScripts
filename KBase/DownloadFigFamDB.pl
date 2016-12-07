use strict;
use Data::Dumper;
use Bio::KBase::DeploymentConfig;
use URI;
use Net::FTP;
use Getopt::Long::Descriptive;

my $default_releasenum = 59;
my $default_release = "Release$default_releasenum";

my($opt, $usage) = describe_options("%c %o [ReleaseNumber ...]",
				    ['help|h', 'Print usage message and exit.'],
				    ['kmer|k=i@', 
				     'Install the given kmer size; repeat the option to install multiple kmer sizes. Defaults to 8mers.', 
				     { default => [ 8 ] }
				    ],
				    [],
				    ["If no release number is supplied the latest production release, $default_releasenum, will be installed."],
				    ["The format of ReleaseNumber is 'ReleaseNN'; e.g. release $default_releasenum is '$default_release'."],
);

print($usage->text), exit  if $opt->help;

my @releases;
if (@ARGV == 0)
{
    push(@releases, $default_release);
}
else
{
    push(@releases, @ARGV);
}

my $cfg = Bio::KBase::DeploymentConfig->new('KmerAnnotationByFigfam',
					    { 
						'kmer-ftp-site' => 'ftp://ftp.theseed.org/FIGfams',
					    });

my $kmer_data = $cfg->setting('kmer-data');
$kmer_data or die "$0: configuration variable kmer-data was not set";
if (! -d $kmer_data)
{
    mkdir $kmer_data || die "$0: cannot mkdir $kmer_data: $!";
}

if (! -d "$kmer_data/ACTIVE")
{
    mkdir "$kmer_data/ACTIVE" || die "$0: cannot mkdir $kmer_data/ACTIVE: $!";
}

my $kmer_tmp = "$kmer_data/tmp";
if (! -d $kmer_tmp)
{
    mkdir $kmer_tmp || die "$0: cannot mkdir $kmer_tmp: $!";
}

    
my $url = $cfg->setting('kmer-ftp-site');
$url or die "$0: configuration variable kmer-ftp-site was not set";
my $uri = URI->new($url);

$uri->scheme eq 'ftp' or die "$0: Only ftp URLs are allowed for kmer-ftp-site";

my $ftp = Net::FTP->new(Host => $uri->host,
			Port => $uri->port);

my $user = $uri->user || "anonymous";
my $hostname = `hostname`;
chomp $hostname;
my $pass = $uri->password || "kbaseuser\@$hostname";

$ftp->login($user, $pass) || die "$0: FTP login failed";

$ftp->cwd($uri->path) || die "$0: FTP cwd " . $uri->path . " failed";

my @files = $ftp->ls();
@files or die "$0: No files found";

my @figfam_files = grep { /^Release\d+\.figfams\.tgz/ } @files;

my %kmers;
my %ffs;
for my $f (@files)
{
    if ($f =~ /^Release(\d+)\.figfams\.tgz/)
    {
	$ffs{$1} = $f;
    }
    elsif ($f =~ /^Release(\d+)\.kmers\.(\d+)\.tgz/)
    {
	$kmers{$1}->{$2} = $f;
    }
}

#print Dumper(\%ffs, \%kmers);

#
# Use the list of requested releases to ensure we have the data we're looking for.
#

my @work;
for my $rel (@releases)
{
    my($n) = $rel =~ /^Release(\d+)/;
    $n or die "$rel is invalid - not of the form ReleaseNN\n";

    my $fams = $ffs{$n};
    $fams or die "Figfams download not available for $rel\n";

    my $kmers = $kmers{$n};
    if (!ref($kmers))
    {
	die "Kmers not available for download for $rel\n";
    }
    my @kwork;
    for my $kmer (@{$opt->kmer})
    {
	if ($kmers->{$kmer})
	{
	    push(@kwork, [$kmer, $kmers->{$kmer}]);
	}
	else
	{
	    die "Kmer size $kmer not available for $rel\n";
	}
    }
    push(@work, [$rel, $n, $fams, \@kwork]);
}
#die Dumper(\@work);

#
# @work is a list of tuples [ReleaseNN, NN, fams-file, $kmer_work]
#
# $kmer_work is a ref to a list of tuples [$kmer-size, $file]
#

for my $ent (@work)
{

    my($rel, $fam, $fam_file, $kmer_work) = @$ent;

    #
    # This code relies on the internal structure of the figfam releases.
    # That's OK.
    #
    if (! -s "$kmer_data/Release$fam/families.2c")
    {
	print "Load base release for $fam\n";
	download_and_unpack($fam_file);
    }
    for my $kent (@$kmer_work)
    {
	my($k, $kmer_file) = @$kent;
	if (! -s "$kmer_data/Release$fam/Merged/$k/table.binary")
	{
	    print "Load k=$k kmers for $fam\n";
	    download_and_unpack($kmer_file);
	}
    }
    if (! -d "$kmer_data/ACTIVE/Release$fam")
    {
	symlink("../Release$fam", "$kmer_data/ACTIVE/Release$fam") or die "symlink ../Release$fam $kmer_data/ACTIVE/Release$fam failed: $!";
    }
}

#
# If no default symlink set, change to the highest release number.
#
if (! -d "$kmer_data/DEFAULT")
{
    opendir(D, $kmer_data) or die "opendir $kmer_data failed: $!";
    my @rels;
    while (my $d = readdir(D))
    {
	if ($d =~ /^Release(\d+)$/)
	{
	    push(@rels, $1);
	}
    }
    closedir(D);
    if (@rels)
    {
	@rels = sort { $b <=> $a } @rels;
	my $max = $rels[0];
	print "Marking release $max as default\n";
	symlink("Release$max", "$kmer_data/DEFAULT") || die "symlink Release$max $kmer_data/DEFAULT failed: $!";
    }
}

sub download_and_unpack
{
    my($file) = @_;
    print "Download $file\n";
    my $file_url = "$url/$file";
    my $rc;
    if (! -s "$kmer_tmp/$file")
    {
	$rc = system("curl", "-L", "-o", "$kmer_tmp/$file", $file_url);
	$rc == 0 or die "Error downloading $file_url to $kmer_tmp/$file";
    }
    print "Download complete; extracting to $kmer_data\n";
    $rc = system("tar", "-C", $kmer_data, "-x", "-z", "-p", "-f", "$kmer_tmp/$file");
    $rc == 0 or die "Error untarring $kmer_tmp/$file";
}