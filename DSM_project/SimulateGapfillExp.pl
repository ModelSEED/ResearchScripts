#!/usr/bin/env perl
########################################################################
# Authors: Christopher Henry, Scott Devoid, Paul Frybarger
# Contact email: chenry@mcs.anl.gov
# Development location: Mathematics and Computer Science Division, Argonne National Lab
########################################################################
use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(save_workspace_object get_workspace_object fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $inputdir = $ARGV[0];

my $values = [
0.0064,
0.0032,
0.0016,
0.0008,
0.0004,
0.0002,
0.0001,
0.00005,
0.000025
];

my $expressionexps = [qw(
S4T1
S101T1
S102T2
S103T1
S104T2
S105T1
S106T1
)];

my $fba = get_fba_client();

my $params = {
          'completeGapfill' => 0,
          'gapFill_workspace' => 'chenry:1431835409789',
          'model' => 'iBsuDSM_isoprene',
          'timePerSolution' => 43200,
          'source_model' => undef,
          'probabilisticAnnotation_workspace' => 'chenry:1431835409789',
          'high_expression_threshold' => '0.00023',
          'high_expression_penalty_factor' => 1,
          'low_expression_penalty_factor' => 1,
          'omega' => 0,
          'singletranspen' => 1,
          'integrate_solution' => 0,
          'nodeltagpen' => 1,
          'biomasstranspen' => 1,
          'workspace' => 'chenry:1431835409789',
          'allowunbalanced' => 0,
          'nostructpen' => 1,
          'directionpen' => 1,
          'gapFill' => undef,
          'drainpen' => 1,
          'num_solutions' => 1,
          'exp_raw_data' => {},
          'expsamplews' => 'chenry:1431835409789',
          'nogprhyp' => 0,
          'unfavorablepen' => 1,
          'transpen' => 1,
          'low_expression_threshold' => '0.00023',
          'wsurl' => 'http://kbase.us/services/ws',
          'solver' => 'cplex',
          'source_model_ws' => 'chenry:1431835409789',
          'formulation' => {
                             'additionalcpds' => [],
                             'nobiomasshyp' => 1,
                             'maximizeObjective' => 1,
                             'bounds' => [],
                             'rxnko' => [],
                             'objfraction' => '0.1',
                             'geneko' => [],
                             'constraints' => [],
                             'uptakelim' => {}
                           },
          'nobiomasshyp' => 1,
          'fastgapfill' => 1,
          'nomediahyp' => 0,
          'gauranteedrxns' => [],
          'alpha' => '0.00007',
          'out_model' => 'expgapfill',
          'expsample' => undef,
          'totalTimeLimit' => 45000,
          'scalefluxes' => 0,
          'blacklistedrxns' => [],
          'nopathwayhyp' => 0,
          'expression_threshold_type' => 'AbsoluteThreshold',
          'model_workspace' => 'chenry:1431835409789',
          'target_reactions' => []
        };

open(my $fhh, ">", $inputdir."/GapfillSummary.txt");
print $fhh "Dataset\tThreshold\tHighexprxn\tLowexprxn\tGapfill\tActive on\tInactive on\tActive off\n";
for (my $j=0; $j < @{$expressionexps}; $j++) {
	$params->{exp_raw_data} = {};
	open(my $fh, "<", $inputdir."/".$expressionexps->[$j].".txt");
	while (my $line = <$fh>) {
		chomp($line);
		my $array = [split(/\t/,$line)];
		$params->{exp_raw_data}->{$array->[0]} = $array->[1];
	}
	close($fh);	
	for (my $i=0; $i < @{$values}; $i++) {
		$params->{low_expression_threshold} = $values->[$i];
		$params->{high_expression_threshold} = $values->[$i];
		$params->{out_model} = "ExpGapfill_".$expressionexps->[$j]."_".$i;
		$fba->gapfill_model($params);
		(my $data,my $prov) = get_workspace_object("chenry:1431835409789/".$params->{out_model}.".gf");
		print $fhh $expressionexps->[$j]."\t".$values->[$i];
		for (my $k=0; $k < 6; $k++) {
			my $array = [split(/:/,$data->{outputfiles}->{gapfillstats}->[$k])];
			print $fhh "\t".$array->[1];
		}
		print $fhh "\n";
		delete $params->{expsample};
	}
	exit;
}
close($fhh);