use Data::Dumper;
use Bio::P3::Workspace::ScriptHelpers;

my $ws = Bio::P3::Workspace::ScriptHelpers::wsClient();
my $username = Bio::P3::Workspace::ScriptHelpers::user();

#1.   Request a list of workspaces shared with a user with limit/offset
my $output = $ws->ls({
	paths => ["/"],
	query => {
		'$or' => [
			{"permissions.".$username => "r"},
			{"permissions.".$username => "w"},
			{"permissions.".$username => "a"}
		]
	}
});
print "1:".Data::Dumper->Dump([$output]);

#2.   Request a list of public workspaces with limit/offset
$output = $ws->ls({
	paths => ["/"],
	query => {
		global_permission => ["r","p"],
	}
});
print "2:".Data::Dumper->Dump([$output]);

#3.   Request list of workspaces owned by a specific user and visible to the requesting user.
$output = $ws->ls({
	paths => ["/reviewer/"]
});
print "3:".Data::Dumper->Dump([$output]);

#4.   Clone an entire workspace (shared, public, or published) to a new workspace owned by the requesting user
$output = $ws->copy({
	objects => [["/reviewer/clone_demo","/".$username."/clone_target"]],
	overwrite => 1,
	recursive => 1
});
print "4:".Data::Dumper->Dump([$output]);