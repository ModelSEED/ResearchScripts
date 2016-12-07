#!/usr/bin/perl -w

use strict;
use YAML;
use Data::Dumper;

my $directory = "/Users/chenry/code/fba_tools/ui/narrative/methods/*";
my $files = [glob($directory)];

my $parameters;
my $methods;
for (my $i=0; $i < @{$files}; $i++) {
	print $files->[$i]."/display.yaml\n";
	my $method;
	if ($files->[$i] =~ /methods\/(.+)/) {
		$method = $1;
	}
	my $data = [YAML::LoadFile($files->[$i]."/display.yaml")];
	foreach my $param (keys(%{$data->[0]->{parameters}})) {
		$parameters->{$param} = $data->[0]->{parameters}->{$param};
		$methods->{$method}->{$param} = 1;
	}
}

foreach my $param (keys(%{$parameters})) {
	chomp($parameters->{$param}->{"ui-name"});
	chomp($parameters->{$param}->{"short-hint"});
	chomp($parameters->{$param}->{placeholder});
	print $param."\t".$parameters->{$param}->{"ui-name"}."\t".$parameters->{$param}->{"short-hint"}."\t".$parameters->{$param}->{placeholder}."\n";
}

foreach my $method (keys(%{$methods})) {
	print $method."\n";
	foreach my $param (keys(%{$methods->{$method}})) {
		print $param."\n";
	}
	print "\n";
}