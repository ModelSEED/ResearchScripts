use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;

my $procs = $ARGV[0];
my $index = $ARGV[1];

my $impl = fba_tools::fba_toolsImpl->new();

if (!defined($ENV{'KB_AUTH_TOKEN'})) {
	require "Bio/KBase/fbaModelServices/ScriptHelpers.pm";
	$ENV{'KB_AUTH_TOKEN'} = Bio::KBase::fbaModelServices::ScriptHelpers::getToken();
}
my $testctx = LocalCallContext->new($ENV{'KB_AUTH_TOKEN'},"chenry",[{'service' => 'fba_tools', 'method' => "unknown", 'method_params' => [{}]}],"unknown");
$fba_tools::fba_toolsServer::CallContext = $testctx;
my $ws = $impl->util_kbase_store()->workspace();
my $models = $ws->list_objects({
	workspaces => ["PubSEEDGramNegModels"],
	type => "KBaseFBA.FBAModel",
});
my $fbas = $ws->list_objects({
	workspaces => ["BiomassSensitivityAnalysis"],
	type => "KBaseFBA.FBA",
});
my $fbahash = {};
for (my $i=0; $i < @{$fbas}; $i++) {
	$fbahash->{$fbas->[$i]->[1]} = 1;
}
my $count = 0;
for (my $i=0; $i < @{$models}; $i++) {
	if ($i % $procs  == $index) {
		#$models->[$i]->[1] = "243276.5.model";
		if (!defined($fbahash->{$models->[$i]->[1].".fba"})) {
			if ($count % 100 == 0) {
				$impl = fba_tools::fba_toolsImpl->new();		}
			eval {
				$impl->gapfill_metabolic_model({
					workspace => "PubSEEDGramNegModels",
					fbamodel_id => $models->[$i]->[1],
					media_id => "Carbon-D-Glucose",
					media_workspace => "KBaseMedia"				
				});
				$impl->run_flux_balance_analysis({
					workspace => "BiomassSensitivityAnalysis",
					fbamodel_id => $models->[$i]->[1],
					fba_output_id => $models->[$i]->[1].".fba",
					fbamodel_workspace => "PubSEEDGramNegModels",
					media_id => "Carbon-D-Glucose",
					media_workspace => "KBaseMedia",
					minimize_flux => 1,
					sensitivity_analysis => 1
				});
			};
			if ($@) {
				my $error = $@;
				print $i.":".$models->[$i]->[1].":".$error."\n";
			} else {
				print $i.":".$models->[$i]->[1].":success\n";
			}
			$count++;
			#exit();
		}
	}
}

{
    package LocalCallContext;
    use strict;
    sub new {
        my($class,$token,$user,$provenance,$method) = @_;
        my $self = {
            token => $token,
            user_id => $user,
            provenance => $provenance,
            method => $method
        };
        return bless $self, $class;
    }
    sub user_id {
        my($self) = @_;
        return $self->{user_id};
    }
    sub token {
        my($self) = @_;
        return $self->{token};
    }
    sub provenance {
        my($self) = @_;
        return $self->{provenance};
    }
    sub method {
        my($self) = @_;
        return $self->{method};
    }
    sub authenticated {
        return 1;
    }
    sub log_debug {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
    sub log_info {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
}