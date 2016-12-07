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

my $cursor = $db->get_collection('workspaces')->find({
	name => "home",
	owner => "nconrad\@patricbrc.org"
});
my $object = $cursor->next;

$cursor = $db->get_collection('objects')->find({
	workspace_uuid => $object->{uuid},
	path => qr/^models\//
});
while ($object = $cursor->next) {
	print $object->{path}."/".$object->{name}."\n";
}

$db->get_collection('objects')->remove({
	workspace_uuid => $object->{uuid},
});