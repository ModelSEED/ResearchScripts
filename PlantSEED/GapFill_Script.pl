#!/usr/bin/env perl
use strict;
use warnings;

my ($ws,$id,$med,$new_id) = ($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
exit if !$ws || !$id;
$med="Hetero" if !$med;
$new_id="TestGF" if !$new_id;

print "Starting: ".scalar(localtime)."\n";

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
#						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
						       'mfatoolkitbin' => "/homes/seaver/Software/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://140.221.85.86:7073/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});

$FBAImpl->_setContext(undef,{auth=>$AToken});
print "Initialized: ".scalar(localtime)."\n";

#######################################################
#Biochemistry and Model
#######################################################

my $modelObj = $FBAImpl->_get_msobject("FBAModel",$ws,$id);
print "Model has ".scalar(@{$modelObj->modelreactions()})." reactions\n";

#######################################################
#MFA
#######################################################

eval {
    $FBAImpl->gapfill_model({workspace=>$ws,model=>$id,out_model=>$new_id,
			     solver=>"CPLEX",completeGapfill=>0,integrate_solution=>0,fastgapfill=>1,
			     formulation=>{formulation=>{media => "Plant".$med."trophicMedia",
							 media_workspace => "seaver:PlantRecModels"},
					   timePerSolution=>87000,
					   totalTimeLimit=>87000}});
};
if ($@) {
    print "Failed at ".scalar(localtime)."\n";
    print "ERROR_MESSAGE".$@."END_ERROR_MESSAGE\n";
}else{
    print "GapFill Finished: ".scalar(localtime)."\n";
}
