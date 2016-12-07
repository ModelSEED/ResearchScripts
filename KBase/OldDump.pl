#!/usr/bin/perl -w

use strict;
use MongoDB::Connection;

$|=1;

my $config = {
	host => "mongodb://wsserv:kbws4ever!\@db3.chicago.kbase.us/workspace_service",
	db_name => "workspace_service",
	auto_connect => 1,
	auto_reconnect => 1
};
my $conn = MongoDB::Connection->new(%$config);
die "Unable to connect: $@" if (!defined($conn));
my $db = $conn->get_database("workspace_service");
my $grid = $db->get_gridfs();
my $cursor = $db->get_collection( 'workspaces' )->find({});
my $workspaces = [];
while (my $obj = $cursor->next) {
	if (!-e $obj->{id}.".tgz") {
		my $objs = $obj->{objects};
        	my $uuids = [];
        	if (@{$uuids} > 0) {
			foreach my $id (keys(%{$objs})) {
        			push(@{$uuids},$objs->{$id});
        		}
        		push(@{$workspaces},[$obj->{id},$uuids]);
		}
	}
}
my $wscount = 0;
my $objcount = 0;
for (my $i=0; $i < @{$workspaces}; $i++) {
        last;
	print $workspaces->[$i]->[0]."\t".$wscount."\t".$objcount."\n";
	$wscount++;
	my $uuids = $workspaces->[$i]->[1];
	system("mkdir ".$workspaces->[$i]->[0]);
	my $objects = [];
	$cursor = $db->get_collection( 'workspaceObjects' )->find({uuid => {'$in' => $uuids}});
	while (my $obj = $cursor->next) {
		push(@{$objects},[$obj->{id},$obj->{chsum},$obj->{type}]);	
	}
	for (my $i=0; $i < @{$objects}; $i++) {
    		$objcount++;
		if ($objects->[$i]->[2] ne "GenomeContigs") {
        		my $file = $grid->find_one({chsum => $objects->[$i]->[1]});
        		if (defined($file)) {
        			my $dataString = $file->slurp();
        			open(my $fh, ">", "NO_WORKSPACE/".$objects->[$i]->[0]);
        			print $fh $dataString;
        			close($fh);
			}
		}
	}
	system("tar -czf ".$workspaces->[$i]->[0].".tgz ".$workspaces->[$i]->[0]);
	system("rm -rf ".$workspaces->[$i]->[0]);
}
$cursor = $db->get_collection( 'workspaceObjects' )->find({workspace => "NO_WORKSPACE",id => ""});
print "NO_WORKSPACE\t".$wscount."\t".$objcount."\n";
system("mkdir NO_WORKSPACE");
my $objects = [];
while (my $obj = $cursor->next) {
	print "Workspace:".$obj->{id}."\n";
	if ($obj->{id} eq "224308.1.anno") {
		$objects = [[$obj->{id},$obj->{chsum},$obj->{type}]];
		last;
	}
	#push(@{$objects},[$obj->{id},$obj->{chsum},$obj->{type}]);
}
for (my $i=0; $i < @{$objects}; $i++) {
    if ($objects->[$i]->[2] ne "GenomeContigs") {
    	print "Printing:".$objects->[$i]->[0]."\n";
    	my $file = $grid->find_one({chsum => $objects->[$i]->[1]});
		if (defined($file)) {
			my $dataString = $file->slurp();
			open(my $fh, ">", "NO_WORKSPACE/".$objects->[$i]->[0]);
			print $fh $dataString;
			close($fh);
		}
    }
}
#system("tar -czf NO_WORKSPACE.tgz NO_WORKSPACE");
#system("rm -rf NO_WORKSPACE");

1;