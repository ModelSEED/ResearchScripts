use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;

my $model = $ARGV[0];#ws/id
my $media = $ARGV[1];#ws/id
my $fba = $ARGV[2];

$model = [split(/\//,$model)];
$media = [split(/\//,$media)];
$fba = [split(/\//,$fba)];

my $impl = fba_tools::fba_toolsImpl->new();
Bio::KBase::utilities::create_context_from_client_config();
$impl->run_flux_balance_analysis({
	workspace => $fba->[0],
	fbamodel_id => $model->[1],
	fba_output_id => $fba->[1],
	fbamodel_workspace => $model->[0],
	media_id => $media->[1],
	media_workspace => $media->[0],
	minimize_flux => 1,
	sensitivity_analysis => 1
});