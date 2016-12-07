use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use File::Temp;
use LWP::Simple;
use Bio::KBase::ObjectAPI::utilities;

my $directory = $ARGV[0];

my $ws = Bio::P3::Workspace::WorkspaceClient->new("https://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');

open (my $fh, "<", $directory."/media");
my $media = "";
while (my $line = <$fh>) {
    $media .= $line;
}
close($fh);
my $mediadata = Bio::KBase::ObjectAPI::utilities::FROMJSON($media);
open ($fh, "<", $directory."/default");
my $biochem = "";
while (my $line = <$fh>) {
    $biochem .= $line;
}
close($fh);
$biochem = Bio::KBase::ObjectAPI::utilities::FROMJSON($biochem);
my $cpdhash = {};
for (my $i=0; $i < @{$biochem->{compounds}}; $i++) {
	$cpdhash->{$biochem->{compounds}->[$i]->{id}} = $biochem->{compounds}->[$i]->{name};
}

foreach my $medianame (keys(%{$mediadata})) {
	my $meta = {
		type => $mediadata->{$medianame}->{type},
		isDefined => $mediadata->{$medianame}->{isDefined},
		source_id => $mediadata->{$medianame}->{source_id},
		isMinimal => $mediadata->{$medianame}->{isMinimal},
		name => $mediadata->{$medianame}->{name},
	};
	my $data = "id\tname\tconcentration\tminflux\tmaxflux\n";
	for (my $i=0; $i < @{$mediadata->{$medianame}->{mediacompounds}}; $i++) {
		if ($mediadata->{$medianame}->{mediacompounds}->[$i]->{compound_ref} =~ m/(cpd\d+)/) {
			my $name = $cpdhash->{$1};
			$data .= $1."\t".
				$name."\t".
				$mediadata->{$medianame}->{mediacompounds}->[$i]->{concentration}."\t".
				$mediadata->{$medianame}->{mediacompounds}->[$i]->{minFlux}."\t".
				$mediadata->{$medianame}->{mediacompounds}->[$i]->{maxFlux}."\n";
		}
	}
	my $list = $ws->create({
		objects => [
			[
				"/chenry/public/modelsupport/media/".$medianame,
				"media",
				$meta,
				$data
			]
		],
		permission => "r",
		createUploadNodes => 0,
		overwrite => 1,
		adminmode => 0,
	});
}