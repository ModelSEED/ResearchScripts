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
my $cursor = $db->get_collection( 'workspaces' )->find({});
my $wscount = 0;
my $objcount = 0;
my $grid = $db->get_gridfs();
while (my $object = $cursor->next) {
        print $object->{id}."\t".$wscount."\t".$objcount."\n";
        $wscount++;
        if (!-e $object->{id}.".tgz") {
	        my $objs = $object->{objects};
	        my $uuids = [];
	        foreach my $id (keys(%{$objs})) {
	                push(@{$uuids},$objs->{$id});
	        }
	        if (@{$uuids} > 0) {
	        	system("mkdir ".$object->{id});
				my $objcursor = $db->get_collection( 'workspaceObjects' )->find({uuid => {'$in' => $uuids}});
				while (my $obj = $objcursor->next) {
					$objcount++;
					if ($obj->{type} ne "GenomeContigs") {
						my $file = $grid->find_one({chsum => $obj->{chsum}});
						if (defined($file)) {
							my $dataString = $file->slurp();
							open(my $fh, ">", $object->{id}."/".$obj->{id});
							print $fh $dataString;
							close($fh);
						}
					}
				}
				system("tar -czf ".$object->{id}.".tgz ".$object->{id});
				system("rm -rf ".$object->{id});
	        }
        }
}
my $objcursor = $db->get_collection( 'workspaceObjects' )->find({workspace => "NO_WORKSPACE"});
print "NO_WORKSPACE\t".$wscount."\t".$objcount."\n";
system("mkdir NO_WORKSPACE");
while (my $obj = $objcursor->next) {
    if ($obj->{type} ne "GenomeContigs") {
    	my $file = $grid->find_one({chsum => $obj->{chsum}});
		if (defined($file)) {
			my $dataString = $file->slurp();
			open(my $fh, ">", "NO_WORKSPACE/".$obj->{id});
			print $fh $dataString;
			close($fh);
		}
    }
}
system("tar -czf NO_WORKSPACE.tgz NO_WORKSPACE");
system("rm -rf NO_WORKSPACE");

1;