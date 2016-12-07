use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use Bio::KBase::AuthToken;

my $ws = Bio::P3::Workspace::WorkspaceClient->new("http://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');
my $wstwo = Bio::P3::Workspace::WorkspaceClient->new("http://p3c.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');

my $wsls = $ws->ls({
	paths => ["/chenry/"],
	recursive => 1,
});
$wsls = $wsls->{"/chenry/"};
for (my $j=0; $j < @{$wsls}; $j++) {
	$wslss = $ws->ls({
		paths => [$wsls->[$j]->[2].$wsls->[$j]->[0]],
		recursive => 1,
	});
	$wslss = $wslss->{$wsls->[$j]->[2].$wsls->[$j]->[0]};
	for (my $k=0; $k < @{$wslss}; $k++) {
		my $output = $ws->get({
			objects => [$wslss->[$k]->[2].$wslss->[$k]->[0]]
		});
		$wstwo->create({
			objects => [[
				$wslss->[$k]->[2].$wslss->[$k]->[0],
				$wslss->[$k]->[1],
				$wslss->[$k]->[7],
				$output->[1],
			]],
			permission => $wslss->[$k]->[10],
			overwrite => 1
		});
	}
}