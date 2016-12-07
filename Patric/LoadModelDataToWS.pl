use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use File::Temp;
use LWP::Simple;
use Bio::KBase::ObjectAPI::utilities;

my $directory = $ARGV[0];

my $ws = Bio::P3::Workspace::WorkspaceClient->new("https://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');

my ($fh1, $classifierFile) = File::Temp::tempfile();
close($fh1);
my $status = LWP::Simple::getstore("http://bioseed.mcs.anl.gov/~chenry/ModelSEED/classifier.txt", $classifierFile);
open (my $fh, "<", $classifierFile);
my $classifierFile = "";
while (my $line = <$fh>) {
    $classifierFile .= $line;
}
close($fh);

open ($fh, "<", $directory."/default");
my $biochem = "";
while (my $line = <$fh>) {
    $biochem .= $line;
}
close($fh);

open ($fh, "<", $directory."/default-mapping");
my $mapping = "";
while (my $line = <$fh>) {
    $mapping .= $line;
}
close($fh);

my $list = $ws->create({
	objects => [
		[
			"/chenry/public/modelsupport/classifiers/gramclassifier.string",
			"string",
			{},
			$classifierFile
		],[
			"/chenry/public/modelsupport/biochemistry/default.biochem",
			"biochemistry",
			{},
			$biochem
		],[
			"/chenry/public/modelsupport/biochemistry/default.mapping",
			"mapping",
			{},
			$mapping
		]
	],
	permission => "r",
	createUploadNodes => 0,
	overwrite => 1,
	adminmode => 0,
});

open ($fh, "<", $directory."/CoreModelTemplate");
my $template = "";
while (my $line = <$fh>) {
    $template .= $line;
}
close($fh);

my $core = Bio::KBase::ObjectAPI::utilities::FROMJSON($template);
$core->{mapping_ref} = $list->[2]->[4];
$core->{biochemistry_ref} = $list->[1]->[4];

open ($fh, "<", $directory."/GramNegModelTemplate");
$template = "";
while (my $line = <$fh>) {
    $template .= $line;
}
close($fh);

my $grampos = Bio::KBase::ObjectAPI::utilities::FROMJSON($template);
$grampos->{mapping_ref} = $list->[2]->[4];
$grampos->{biochemistry_ref} = $list->[1]->[4];

open ($fh, "<", $directory."/GramPosModelTemplate");
$template = "";
while (my $line = <$fh>) {
    $template .= $line;
}
close($fh);

my $gramneg = Bio::KBase::ObjectAPI::utilities::FROMJSON($template);
$gramneg->{mapping_ref} = $list->[2]->[4];
$gramneg->{biochemistry_ref} = $list->[1]->[4];

$list = $ws->create({
	objects => [
		[
			"/chenry/public/modelsupport/templates/Core.template",
			"modeltemplate",
			{},
			Bio::KBase::ObjectAPI::utilities::TOJSON($core)
		],[
			"/chenry/public/modelsupport/templates/GramPositive.template",
			"modeltemplate",
			{},
			Bio::KBase::ObjectAPI::utilities::TOJSON($grampos)
		],[
			"/chenry/public/modelsupport/templates/GramNegative.template",
			"modeltemplate",
			{},
			Bio::KBase::ObjectAPI::utilities::TOJSON($gramneg)
		]
	],
	permission => "r",
	createUploadNodes => 0,
	overwrite => 1,
	adminmode => 0,
});