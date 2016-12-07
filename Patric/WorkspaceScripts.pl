use Data::Dumper;
use Bio::P3::Workspace::WorkspaceClient;
use Bio::P3::Workspace::WorkspaceImpl;
use Bio::KBase::AuthToken;

my $ws = Bio::P3::Workspace::WorkspaceClient->new("https://p3.theseed.org/services/Workspace",user_id => 'chenry', password => 'hello824');
$Bio::P3::Workspace::Service::CallContext = Bio::P3::Workspace::ServiceContext->new("un=chenry|tokenid=03B0C858-7A70-11E4-9DE6-FDA042A49C03|expiry=1449094224|client_id=chenry|token_type=Bearer|SigningSubject=http://rast.nmpdr.org/goauth/keys/E087E220-F8B1-11E3-9175-BD9D42A49C03|sig=085255b952c8db3ddd7e051ac4a729f719f22e531ddbc0a3edd86a895da851faa93249a7347c75324dc025b977e9ac7c4e02fb4c966ec6003ecf90d3148e35160265dbcdd235658deeed0ec4e0c030efee923fda1a55e8cc6f116bcd632fa6a576d7bf4a794554d2d914b54856e1e7ac2b071f81a8841d142123095f6af957cc","test","chenry");
$ENV{KB_DEPLOYMENT_CONFIG}="/Users/chenry/code/Workspace/configs/test.cfg";
my $wslocal = Bio::P3::Workspace::WorkspaceImpl->new();

my $wsls = $ws->ls({
	paths => ["/chenry/models/"],
	recursive => 1,
	adminmode => 1
});
$wsls = $wsls->{"/chenry/models/"};
$wsls = [sort { $a->[2] cmp $b->[2] } @{$wsls}];
for (my $j=0; $j < @{$wsls}; $j++) {
	print $wsls->[$j]->[2].$wsls->[$j]->[0]."\t".$wsls->[$j]->[1]."\n";
}

#$wslocal->delete({
#	objects => ["/chenry/models/.TestModel"],
#	deleteDirectories => 1,
#	force => 1
#});

package Bio::P3::Workspace::ServiceContext;

use strict;

sub new {
    my($class,$token,$method,$user) = @_;
    my $self = {
        token => $token,
        method => $method,
        user_id => $user
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
sub method {
	my($self) = @_;
	return $self->{method};
}