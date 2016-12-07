#!/usr/bin/perl -w

use strict;
use warnings;
use DBI;
use Config::Simple;
use Bio::KBase::workspaceService::Client;
use Bio::KBase::fbaModelServices::Client;
$|=1;

my $config = $ARGV[0];
if (!defined($config)) {
	print STDERR "No config file provided!\n";
	exit(-1);
}
if (!-e $config) {
	print STDERR "Config file ".$config." not found!\n";
	exit(-1);
}
#Params: kbclientconfig.wsurl, kbclientconfig.fbaurl, kbclientconfig.auth
my $c = Config::Simple->new();
$c->read($config);

my $wserv = Bio::KBase::workspaceService::Client->new($c->param("kbclientconfig.wsurl"));

my $jobs = $wserv->get_jobs({
	type => "FBA",
	status => "error",
	auth => $c->param("kbclientconfig.auth")
});

print STDERR @{$jobs}." jobs with errors!\n";
print "ID\tGenome\tOwner\tError type\tDetails\n";
for (my $i=0; $i < @{$jobs}; $i++) {
	my $job = $jobs->[$i];
	my $resubmit = 0;
	my $delete = 0;
	if ($job->{jobdata}->{error} =~ m/Can.t\suse\sstring\s\(.job\.\d+.\)\sas\sa\sHASH\sref\swhile/) {
		print "Local job bug!\n";
		#$resubmit = 1;
	} elsif ($job->{jobdata}->{error} =~ m/Gapfilling\scompleted.\sbut\sno\svalid\ssolutions\sfound\safter\s24\shours./) {
		print "No gapfilling solutions!\n";
	} else {
		print STDERR $job->{id}.":\n".$job->{jobdata}->{error}."\n\n";
		#$resubmit = 1;
	}
	if ($resubmit == 1) {
		eval {
			$wserv->set_job_status({
				auth => $c->param("kbclientconfig.auth"),
				jobid => $job->{id},
				currentStatus => "error",
				status => "queued",
			});
		};
	} elsif ($delete == 1) {
		eval {
			$wserv->set_job_status({
				auth => $c->param("kbclientconfig.auth"),
				jobid => $job->{id},
				currentStatus => "error",
				status => "delete",
			});
		};
	}	
}	

1;
