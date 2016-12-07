#!/usr/bin/perl

use strict;
use REST::Client;
use LWP::UserAgent;
use JSON::XS;
use Data::Dumper;
use HTTP::Request::Common;

my $url = "http://p3.theseed.org/services/shock_api";
my $ua = LWP::UserAgent->new();
#my $rest = REST::Client->new(useragent => $ua);

my $auth_token = "un=chenry|tokenid=03B0C858-7A70-11E4-9DE6-FDA042A49C03|expiry=1449094224|client_id=chenry|token_type=Bearer|SigningSubject=http://rast.nmpdr.org/goauth/keys/E087E220-F8B1-11E3-9175-BD9D42A49C03|sig=085255b952c8db3ddd7e051ac4a729f719f22e531ddbc0a3edd86a895da851faa93249a7347c75324dc025b977e9ac7c4e02fb4c966ec6003ecf90d3148e35160265dbcdd235658deeed0ec4e0c030efee923fda1a55e8cc6f116bcd632fa6a576d7bf4a794554d2d914b54856e1e7ac2b071f81a8841d142123095f6af957cc";
#my @auth_header = ();
#$rest->addHeader(@auth_header);

#Downloading the file
my $rest = REST::Client->new(useragent => $ua);
$rest->addHeader(Authorization => "OAuth $auth_token");
my $res = $rest->GET($url . "/node/fdce4639-0aef-4340-9b20-aa9572f0c15b?download");
if ($rest->responseCode != 200){
	die "get_file failed: " . $rest->responseContent();
}
print $rest->responseContent();

#Uploading a file
#my $req = HTTP::Request::Common::POST($url."/node/fdce4639-0aef-4340-9b20-aa9572f0c15b",Authorization => "OAuth $auth_token",Content_Type => 'multipart/form-data',Content => [upload => ["/Users/chenry/workspace/example.json"]]);
#$req->method('PUT');
#my $res = $ua->request($req);
       
#Ceating the empty node
#my $res = $ua->post($url."/node",Authorization => "OAuth $auth_token");

#Adding the user to the list of ACLs
#my $res = $ua->put($url."/node/fdce4639-0aef-4340-9b20-aa9572f0c15b/acl/all?users=reviewer",Authorization => "OAuth $auth_token");

my $json = JSON::XS->new;
my $ret = $json->decode($res->content);
print Dumper($ret);

1;