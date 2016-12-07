#!/usr/bin/perl -w

use strict;
use Bio::KBase::fbaModelServices::ScriptHelpers qw(get_workspace_object);
use Bio::KBase::workspace::ScriptHelpers qw(workspaceURL get_ws_client workspace parseObjectMeta parseWorkspaceMeta);
use JSON::XS;
use Digest::MD5;
use DateTime;
use Data::Dumper;
$|=1;

my $genehash;
open(my $fh, "<", "/Users/chenry/workspace/KBase Project/MR1.tab");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if ($array->[4] =~ m/(SO_\d+)/) {
		$genehash->{$1} = $array->[3];
	}
}
close($fh);

my $wsClient = get_ws_client();
(my $data,my $info) = get_workspace_object("jplfaria:1464221538068/Shew_GenBank");
for (my $i=0; $i < @{$data->{features}}; $i++) {
	if (defined($genehash->{$data->{features}->[$i]->{id}})) {
		$data->{features}->[$i]->{function} = $genehash->{$data->{features}->[$i]->{id}};
		print $data->{features}->[$i]->{function}."\n";
	}
}

$info = $wsClient->save_objects({
	'id' => 7537,
	'objects' => [{
	    type => 'KBaseGenomes.Genome',
	    data => $data,
	    name => "Shew_Uniprot",
	    'meta' => {},
	    'provenance' => []
	}]
});

1;