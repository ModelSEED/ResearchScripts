use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
use Bio::KBase::ObjectAPI::utilities;

my $directory = $ARGV[0];
my $list = get_ws_objects_list("KBaseMedia");
my $mediadata = {};
for (my $i=0; $i < @{$list}; $i++) {
	my ($data,$info) = get_workspace_object("KBaseMedia/".$list->[$i]->[1]);
	$mediadata->{$list->[$i]->[1]} = $data;
}
open ( my $fh, ">", $directory."/media");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($mediadata);
close($fh);
