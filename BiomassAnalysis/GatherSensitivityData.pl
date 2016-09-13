use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;

local $| = 1;
my $impl = fba_tools::fba_toolsImpl->new();

if (!defined($ENV{'KB_AUTH_TOKEN'})) {
	require "Bio/KBase/fbaModelServices/ScriptHelpers.pm";
	$ENV{'KB_AUTH_TOKEN'} = Bio::KBase::fbaModelServices::ScriptHelpers::getToken();
}
my $testctx = LocalCallContext->new($ENV{'KB_AUTH_TOKEN'},"chenry",[{'service' => 'fba_tools', 'method' => "unknown", 'method_params' => [{}]}],"unknown");
$fba_tools::fba_toolsServer::CallContext = $testctx;
my $ws = $impl->util_kbase_store()->workspace();
my $fbas = $ws->list_objects({
	workspaces => ["BiomassSensitivityAnalysis"],
	type => "KBaseFBA.FBA",
});
my $fbahash = {};
my $modelhash = {};
my $biohisto = {};
my $biogfhisto = {};
my $modelcount = 0;
print "FBA count:".@{$fbas}."\n";
my $numbio = 0;
my $numcoupled = 0;
for (my $i=0; $i < @{$fbas}; $i++) {
#for (my $i=0; $i < 100; $i++) {
	my $rxnhash = {};
	my $modelid = substr($fbas->[$i]->[1],0,length($fbas->[$i]->[1])-4);
	my $output = $ws->get_objects([{
		workspace => "BiomassSensitivityAnalysis",
		name => $fbas->[$i]->[1]
	},{
		workspace => "PubSEEDGramNegModels",
		name => $modelid
	}]);
	$modelhash->{$fbas->[$i]->[1]} = {
		numbio => 0,
		numcoupled => 0,
		biohash => {},
		biogfhash => {}
	};
	my $rxns = $output->[0]->{data}->{FBAReactionVariables};
	for (my $j=0; $j < @{$output->[1]->{data}->{modelreactions}}; $j++) {
		$rxnhash->{$output->[1]->{data}->{modelreactions}->[$j]->{id}} = $output->[1]->{data}->{modelreactions}->[$j];
	}
	for (my $j=0; $j < @{$rxns}; $j++) {
		my @array = split(/\//,$rxns->[$j]->{modelreaction_ref});
		my $id = pop(@array);
		if (defined($rxns->[$j]->{coupled_reactions}) && @{$rxns->[$j]->{coupled_reactions}} > 0) {
			$modelhash->{$fbas->[$i]->[1]}->{numcoupled}++;
		}
		if (defined($rxns->[$j]->{biomass_dependencies}) && @{$rxns->[$j]->{biomass_dependencies}} > 0) {
			$modelhash->{$fbas->[$i]->[1]}->{numbio}++;
			for (my $k=0; $k < @{$rxns->[$j]->{biomass_dependencies}}; $k++) {
				if (!defined($modelhash->{$fbas->[$i]->[1]}->{biohash}->{$rxns->[$j]->{biomass_dependencies}->[$k]->[1]})) {
					$modelhash->{$fbas->[$i]->[1]}->{biohash}->{$rxns->[$j]->{biomass_dependencies}->[$k]->[1]} = 0;
					if (defined($rxnhash->{$id}) && keys(%{$rxnhash->{$id}->{gapfill_data}}) > 0) {
						$modelhash->{$fbas->[$i]->[1]}->{biogfhash}->{$rxns->[$j]->{biomass_dependencies}->[$k]->[1]} = 0;
					}
				}
				$modelhash->{$fbas->[$i]->[1]}->{biohash}->{$rxns->[$j]->{biomass_dependencies}->[$k]->[1]}++;
				if (defined($rxnhash->{$id}) && keys(%{$rxnhash->{$id}->{gapfill_data}}) > 0) {
					$modelhash->{$fbas->[$i]->[1]}->{biogfhash}->{$rxns->[$j]->{biomass_dependencies}->[$k]->[1]}++;
				}
			}
		}
	}
	foreach my $key (keys(%{$modelhash->{$fbas->[$i]->[1]}->{biohash}})) {
		if (!defined($biohisto->{$key})) {
			$biohisto->{$key}->[0] = $modelcount;
		}
		$biohisto->{$key}->[$modelhash->{$fbas->[$i]->[1]}->{biohash}->{$key}]++;
	}
	foreach my $key (keys(%{$biohisto})) {
		if (!defined($modelhash->{$fbas->[$i]->[1]}->{biohash}->{$key})) {
			$biohisto->{$key}->[0]++;
		}
	}
	foreach my $key (keys(%{$modelhash->{$fbas->[$i]->[1]}->{biogfhash}})) {
		if (!defined($biogfhisto->{$key})) {
			$biogfhisto->{$key}->[0] = $modelcount;
		}
		$biogfhisto->{$key}->[$modelhash->{$fbas->[$i]->[1]}->{biogfhash}->{$key}]++;
	}
	foreach my $key (keys(%{$biogfhisto})) {
		if (!defined($modelhash->{$fbas->[$i]->[1]}->{biogfhash}->{$key})) {
			$biogfhisto->{$key}->[0]++;
		}
	}
	if ($modelhash->{$fbas->[$i]->[1]}->{numbio} > 0) {
		$numbio++;
	}
	if ($modelhash->{$fbas->[$i]->[1]}->{numcoupled} > 0) {
		$numcoupled++;
	}
	print $i."\t".$fbas->[$i]->[1]."\t".$modelhash->{$fbas->[$i]->[1]}->{numbio}."\t".$modelhash->{$fbas->[$i]->[1]}->{numcoupled}."\n";
	$modelcount++;
}

print "\n";

foreach my $key (keys(%{$biohisto})) {
	print $key;
	for (my $i=0; $i < 1000; $i++) {
		print "\t".$biohisto->{$key}->[$i];
	}
	print "\n";
}

print "\n";

foreach my $key (keys(%{$biogfhisto})) {
	print $key;
	for (my $i=0; $i < 1000; $i++) {
		print "\t".$biogfhisto->{$key}->[$i];
	}
	print "\n";
}

print "Model";
foreach my $key (keys(%{$biohisto})) {
	print "\t".$key."_gf";
}
foreach my $key (keys(%{$biohisto})) {
	print "\t".$key;
}
print "\n";
foreach my $model (keys(%{$modelhash})) {
	print $model;
	foreach my $key (keys(%{$biohisto})) {
		print "\t";
		if (defined($modelhash->{$model}->{biogfhash}->{$key})) {
			print $modelhash->{$model}->{biogfhash}->{$key};
		} else {
			print 0;
		}
	}
	foreach my $key (keys(%{$biohisto})) {
		print "\t";
		if (defined($modelhash->{$model}->{biohash}->{$key})) {
			print $modelhash->{$model}->{biohash}->{$key};
		} else {
			print 0;
		}
	}	
	print "\n";
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