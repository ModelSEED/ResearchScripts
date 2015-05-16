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
	name => "models",
	owner => "chenry"
});
my $object = $cursor->next;

$db->get_collection('objects')->remove({
	workspace_uuid => $object->{uuid},
});