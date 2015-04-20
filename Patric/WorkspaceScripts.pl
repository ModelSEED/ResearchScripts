use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use Bio::KBase::AuthToken;

my $ws = Bio::P3::Workspace::WorkspaceClient->new("https://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');
my $output = $ws->ls({
	paths => ["/dmachi\@patricbrc.org/home/"],
	adminmode => 1
});
print Data::Dumper->Dump([$output])."\n";