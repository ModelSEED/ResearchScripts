#!/usr/bin/perl -w

use strict;
use Bio::KBase::workspaceService::Impl;
use Bio::KBase::fbaModelServices::Impl;

my $fba = Bio::KBase::fbaModelServices::Impl->new({
	"fba-url" => "",
	"probanno-url" => "",
	"mssserver-url" => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
	accounttype => "kbase",
	defaultJobState => "queued",
	"idserver-url" => "http://bio-data-1.mcs.anl.gov/services/idserver"
});
my $wserv = Bio::KBase::workspaceService::Impl->new({
	accounttype => "kbase",
#	"idserver-url" => "http://bio-data-1.mcs.anl.gov/services/idserver"
});
for (my $i=94; $i < 14909; $i++) {
	my $id = $wserv->_get_new_id("job.");
	print "Job ID:".$id."\n";
}
1;
