#!/usr/bin/perl -w

use strict;
use Config::Simple;
use Bio::KBase::workspaceService::Client;
use JSON::XS;

my $config = $ARGV[0];
my $directory = $ARGV[1];
if (!defined($config)) {
	print STDERR "No config file provided!\n";
	exit(-1);
}
if (!-e $config) {
	print STDERR "Config file ".$config." not found!\n";
	exit(-1);
}

open( my $fh, "<", $directory."workspaces.list");
my $workspaces;
{
    local $/;
    my $str = <$fh>;
    $workspaces = decode_json $str;
}
close($fh);

my $c = Config::Simple->new();
$c->read($config);

my $wserv = Bio::KBase::workspaceService::Client->new($c->param("kbclientconfig.wsurl"));

for (my $i=0; $i < @{$workspaces}; $i++) {
	my $list = $wserv->list_workspace_objects({
		workspace => $workspaces->[$i]->{id},
		auth => $c->param("kbclientconfig.auth")
	});
	for (my $j=0; $j < @{$list}; $j++) {
		print $list->[$j]->[1]."\t".$workspaces->[$i]->{id}."\t".$list->[$j]->[0]."\t".$list->[$j]->[6]."\t".$list->[$j]->[4]."\t".$list->[$j]->[2]."\n";
	}
}

1;
