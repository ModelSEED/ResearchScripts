use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use fba_tools::fba_toolsImpl;

local $| = 1;
if (!defined($ENV{'KB_AUTH_TOKEN'})) {
	require "Bio/KBase/fbaModelServices/ScriptHelpers.pm";
	$ENV{'KB_AUTH_TOKEN'} = Bio::KBase::fbaModelServices::ScriptHelpers::getToken();
}
my $tester = LocalTester->new($ENV{'KB_AUTH_TOKEN'},$ENV{'KB_DEPLOYMENT_CONFIG'},"chenry:1454960620516");

my $procs = $ARGV[0];
my $index = $ARGV[1];

my $ws = $tester->wsclient();
my $genomes = $ws->list_objects({
	workspaces => ["coremodels"],
	type => "KBaseGenomes.Genome",
});

my $count = 0;
for (my $i=0; $i < @{$genomes}; $i++) {
	if ($i % $procs  == $index) {
		if ($count % 100 == 0) {
			$tester = LocalTester->new($ENV{'KB_AUTH_TOKEN'},$ENV{'KB_DEPLOYMENT_CONFIG'},"chenry:1454960620516");
		}
		eval {
			$tester->impl()->build_metabolic_model({
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
	package LocalTester;
	use strict;
	use Test::More;
    sub new {
        my ($class,$token,$configfile,$wsname) = @_;
        my $c = Config::Simple->new($configfile);
        my $config = $c->get_block('fba_tools');
        my $object = fba_tools::fba_toolsImpl->new();
        my $self = {
            token => $token,
            config_file => $configfile,
            config => $config,
            user_id => undef,
            ws_client => undef,
            ws_name => $wsname,
            obj => $object,
            testcount => 0,
            completetestcount => 0,
            dumpoutput => 0,
            testoutput => {},
            showerrors => 1
        };
        my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
        $self->{user_id} = $auth_token->user_id();
        print "Workspace:".$self->{config}->{"workspace-url"}."\n";
        $self->{ws_client} = new Bio::KBase::workspace::Client($self->{config}->{"workspace-url"},token => $self->{token});
        if (!defined($self->{ws_name})) {
	        my $suffix = int(time * 1000);
	        $self->{ws_name} = 'test_fba_tools_' . $suffix;
	        $self->{ws_client}->create_workspace({workspace => $self->{ws_name}});
	    }
	    print "Test output saved to:".$self->{ws_name}."\n";
        return bless $self, $class;
    }
    sub impl {
    	my ($self) = @_;
    	return $self->{obj};
    }
    sub wsclient {
    	my ($self) = @_;
    	return $self->{ws_client};
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