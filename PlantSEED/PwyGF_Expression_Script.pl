#!/usr/bin/env perl
use strict;
use warnings;
my @temp=();

#######################################################
#Initialization
#######################################################

use Bio::KBase::fbaModelServices::ScriptHelpers qw( getToken );
my $AToken = getToken();

use Bio::KBase::fbaModelServices::Impl;
my $FBAImpl = Bio::KBase::fbaModelServices::Impl->new({'fbajobcache' => "/homes/seaver/Projects/KBase_Scripts/FBA_Scripts/JobsCache",
						       'jobserver-url' => "http://kbase.us/services/workspace",
						       'fbajobdir' => "/tmp/fbajobs",
						       'mfatoolkitbin' => "/vol/model-prod/kbase/MFAToolkit/bin/mfatoolkit",
						       'probanno-url' => "http://140.221.85.86:7073/",
						       'mssserver-url' => "http://bio-data-1.mcs.anl.gov/services/ms_fba",
						       'accounttype' => "kbase",
						       'workspace-url' => "http://kbase.us/services/ws",
						       'defaultJobState' => "queued",
						       'gaserver-url' => "http://kbase.us/services/genome_annotation",
						       'idserver-url' => "http://kbase.us/services/idserver"});
$FBAImpl->_setContext(undef,{auth=>$AToken});

my $modelObj = $FBAImpl->_get_msobject("FBAModel","Maize_Models","Evidenced-Maize-Model");

my $WSClient = $FBAImpl->_workspaceServices();
my $sample = $WSClient->get_object({id=>"Maize_qTeller_Mature_leaf",type=>"KBaseExpression.ExpressionSample",workspace=>"seaver:MaizeTissues"});

#print join("\n", keys %{$sample->{"data"}->{"expression_levels"}}),"\n";
#__END__

my $exp_scores = {};
my $no_reaction = {};
my $no_prot = {};
my $no_subunit = {};
my $miss_protein = {};
my $no_ftr = {};
my $yes_ftr = {};
my $unknown = {};
foreach my $mdlrxn (@{$modelObj->modelreactions()}) {
    $exp_scores->{$mdlrxn->id()} = 0;
    $miss_protein->{$mdlrxn->id()} = 1 if (@{$mdlrxn->modelReactionProteins()} == 0);

    $unknown->{$mdlrxn->id()} = 1 if (@{$mdlrxn->modelReactionProteins()} == 0);

    # Maximal gene expression for a reaction
    my $rxn_score = POSIX::FLT_MIN;
    foreach my $prot (@{$mdlrxn->modelReactionProteins()}) {
	
	# Minimal gene expression for a complex
	my $prot_score = POSIX::FLT_MAX;
	foreach my $subunit (@{$prot->modelReactionProteinSubunits()}) {

	    if (@{$subunit->features()} == 0) {
		next; # Not last, since there may be scores for other subunits
	    }
	    
            # Maximal gene expression for a subunit
	    my $subunit_score = POSIX::FLT_MIN;
	    foreach my $feature (@{$subunit->features()}) {

#		print $feature->id(),"\n";
		my $ftr_id = $feature->id();
#		$ftr_id =~ s/\.\d+\.CDS$//;

		if (!exists $sample->{"data"}->{"expression_levels"}->{$feature->id()}) {
		    $no_ftr->{$feature->id()}=1;

                    # exp level is Unknown, so we can't compute a score for the subunit
		    # $subunit_score = POSIX::FLT_MAX;
		    # last;

		}else{

		    my $ftr_score = $sample->{"data"}->{"expression_levels"}->{$feature->id};
		    $yes_ftr->{$feature->id()}=1;

		    $subunit_score = ($subunit_score > $ftr_score) ? $subunit_score : $ftr_score;

		}
	    }
	    
	    if($subunit_score == POSIX::FLT_MIN){
		$no_subunit->{$mdlrxn->id()}=1;
	    }

	    $prot_score = ($prot_score < $subunit_score) ? $prot_score : $subunit_score;
	}

	if ($prot_score == POSIX::FLT_MAX) {
	    # there was insufficient information to compute a protein score
	    $no_prot->{$mdlrxn->id()} = 1;
	}

	$rxn_score = ($rxn_score > $prot_score) ? $rxn_score : $prot_score;
    
    }

    if ($rxn_score == POSIX::FLT_MIN) {
	# there was insufficient information to compute a protein score
	$no_reaction->{$mdlrxn->id()} = 1;
    }

    $exp_scores->{$mdlrxn->id()} = $rxn_score;

}

#print scalar(keys %$no_subunit),"\t",scalar(keys %$no_prot),"\t",scalar(keys %$no_reaction),"\t",scalar(keys %$exp_scores),"\t",scalar(@{$modelObj->modelreactions()}),"\n";
#__END__
#print scalar(keys %$exp_scores),"\t",scalar(keys %$no_score),"\t",scalar(keys %$no_protein),"\n";
#print scalar(keys %$no_ftr),"\t",scalar(keys %$yes_ftr),"\n";
#print join("\n", map { $exp_scores->{$_} } keys %$exp_scores),"\n";

my $final_exp_scores = {};
my $min_exp_score = (sort {$a <=> $b} grep { $_ != POSIX::FLT_MIN && $_ != POSIX::FLT_MAX } map { $exp_scores->{$_} } keys %$exp_scores)[0];
my $max_exp_score = (sort {$b <=> $a} grep { $_ != POSIX::FLT_MIN && $_ != POSIX::FLT_MAX } map { $exp_scores->{$_} } keys %$exp_scores)[0];

foreach my $rxn_id (keys %$exp_scores) {
    if (exists $unknown->{$rxn_id}) {
#	$final_exp_scores->{$rxn_id} = 1; # i.e., 100% flux bounds when expression score is not available
    } else {
	if($exp_scores->{$rxn_id} eq POSIX::FLT_MIN){
	    $final_exp_scores->{$rxn_id}=0;
	}elsif($exp_scores->{$rxn_id} eq POSIX::FLT_MAX){
	    $final_exp_scores->{$rxn_id}=1;
	}else{
	    # Then normalize the picked sample's expression score by the max
	    $final_exp_scores->{$rxn_id} = sprintf("%.4f",($exp_scores->{$rxn_id} / $max_exp_score));
	}
    }
}

#print join("\n", map {$final_exp_scores->{$_}} keys %$final_exp_scores),"\n";

open(FH, "< ".$ENV{SEAVER_PROJECT}."KBase_Scripts/Pathway_GapFill/Tissue_Pwy_GF/InactiveModelReactions.txt");
open(HI, "> ".$ENV{SEAVER_PROJECT}."KBase_Scripts/Pathway_GapFill/High_Expression.txt");
open(LO, "> ".$ENV{SEAVER_PROJECT}."KBase_Scripts/Pathway_GapFill/Low_Expression.txt");

while(<FH>){
    chomp;
    @temp=split(/\t/,$_,-1);
    if($temp[0] eq "bio1"){
	print HI "bio1	1\n";
	next;
    }

    if(exists($final_exp_scores->{$temp[0]})){
	if($final_exp_scores->{$temp[0]} < 0.2){
	    print LO $modelObj->getObject("modelreactions",$temp[0])->direction(),"\t",$temp[0],"\t0.1\n";
	}else{
	    print HI $temp[0],"\t",$final_exp_scores->{$temp[0]},"\n";
	}
    }else{
#	print HI $temp[0],"\t0.1\n";
	print LO $modelObj->getObject("modelreactions",$temp[0])->direction(),"\t",$temp[0],"\t0.1\n";
    }
}
