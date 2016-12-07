use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use File::Temp;
use LWP::Simple;
use LWP::UserAgent;
use Bio::KBase::ObjectAPI::utilities;
use Bio::P3::Workspace::ScriptHelpers;

my $directory = $ARGV[0];
open (my $fh, "<", $directory."/PubSEED.txt");
my $pubseedhash = {};
while (my $line = <$fh>) {
   chomp($line);
   my $array = [split(/\t/,$line)];
   $pubseedhash->{$array->[0]} = $array->[1];
}
close($fh);

open (my $fh, "<", $directory."/ModelList.txt");
my $modelhash = {};
while (my $line = <$fh>) {
   chomp($line);
   my $array = [split(/\t/,$line)];
   $modelhash->{$array->[0]} = $array;
}
close($fh);

my $ws = Bio::P3::Workspace::WorkspaceClient->new("https://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');

foreach my $mdl (keys(%{$modelhash})) {
	my $res;
	eval {
		$res = Bio::P3::Workspace::ScriptHelpers::wscall("get",{
			objects => ["/".$modelhash->{$mdl}->[2]."/modelseed/".$mdl],
			adminmode => 1
		});
	};
	if ($@) {
		print "Missing:".$mdl."\n";
	} else {
		print "Processing:".$mdl."\n";
		if (defined($pubseedhash->{$modelhash->{$mdl}->[1]})) {
			$res->[0]->[0]->[7]->{genome_source} = "PubSEED";
		} else {
			$res->[0]->[0]->[7]->{genome_source} = "RAST";
		}
		$res = Bio::P3::Workspace::ScriptHelpers::wscall("update_metadata",{
			objects => [[
				"/".$modelhash->{$mdl}->[2]."/modelseed/".$mdl,
				$res->[0]->[0]->[7],
				"modelfolder"
			]],
			adminmode => 1
		});
	}
}