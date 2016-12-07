#!/usr/bin/perl

use strict;
use warnings;

open( my $fh, "<", $ARGV[0]);
my $hash;
my $hashgenomes;
while (my $line = <$fh>) {
	if ($line =~ m/The\sJSONRPC\sserver\sinvocation\sof\sthe\smethod\s.getRastGenomeData.\sfailed\swith\sthe\sfollowing\serror./) {
		$line = <$fh>;
		if ($line =~ m/ERROR.(.+).ERROR/) {
			my $item = $1;
			if (!defined($hash->{$item})) {
				$hash->{$item} = 0;
			}
			$hash->{$item}++;
			my $continue = 1;
			while ($continue == 1) {
				$line = <$fh>;
				if ($line =~ m/NoSuchObjectException.\sNo\sobject\swith\sname\s(.+)\sexists\sin\sworkspace\s(\d+)/) {
					$hashgenomes->{$item}->{$2."/".$1} = 1;
					$continue = 0;
				}
			}
		}
	}	
}
close($fh);

foreach my $key (keys(%{$hash})) {
	print $key."\t".$hash->{$key}."\t".join(";",keys(%{$hashgenomes->{$key}}))."\n";
}