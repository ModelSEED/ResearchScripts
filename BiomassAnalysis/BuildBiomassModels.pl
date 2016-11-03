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
my $genomes = $ws->list_objects({
	workspaces => ["coremodels"],
	type => "KBaseGenomes.Genome",
});

my $count = 0;
for (my $i=0; $i < @{$genomes}; $i++) {
	if ($i % $procs  == $index) {
		if ($count % 100 == 0) {
			$impl = fba_tools::fba_toolsImpl->new();		}
		eval {
			$impl->build_metabolic_model({
				workspace => "PubSEEDGramNegModels",
				genome_id => $genomes->[$i]->[1],
				fbamodel_output_id => $genomes->[$i]->[1].".model",
				media_id => "Carbon-D-Glucose",
				template_id => "gramneg",
				genome_workspace => "coremodels",
				media_workspace => "KBaseMedia",
				gapfill_model => 1
			});
		};
		if ($@) {
			my $error = $@;
			print $i.":".$genomes->[$i]->[1].":".$error;
		} else {
			print $i.":".$genomes->[$i]->[1].":success";
		}
		$count++;
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