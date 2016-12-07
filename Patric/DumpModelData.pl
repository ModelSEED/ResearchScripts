use Data::Dumper;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );
use Bio::KBase::ObjectAPI::utilities;

my $directory = $ARGV[0];

my $kbasews = get_ws_client();
my ($data,$info) = get_workspace_object("kbase/default");
my $biochem = $data;
($data,$info) = get_workspace_object("kbase/default-mapping");
my $mapping = $data;
($data,$info) = get_workspace_object("KBaseTemplateModels/GramPosModelTemplate");
my $gpostemp = $data;
($data,$info) = get_workspace_object("KBaseTemplateModels/GramNegModelTemplate");
my $gnegtemp = $data;
($data,$info) = get_workspace_object("KBaseTemplateModels/CoreModelTemplate");
my $coretemp = $data;

open ( my $fh, ">", $directory."/default");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($biochem);
close($fh);
open ( $fh, ">", $directory."/default-mapping");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($mapping);
close($fh);
open ( $fh, ">", $directory."/GramPosModelTemplate");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($gpostemp);
close($fh);
open ( $fh, ">", $directory."/GramNegModelTemplate");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($gnegtemp);
close($fh);
open ( $fh, ">", $directory."/CoreModelTemplate");
print $fh Bio::KBase::ObjectAPI::utilities::TOJSON($coretemp);
close($fh);