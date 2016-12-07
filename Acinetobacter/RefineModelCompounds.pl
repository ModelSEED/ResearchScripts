#!/usr/bin/perl

$|=1;

my $filename = $ARGV[0];

my $linearray = [];
open(my $fh, "<", $filename);
while (my $line = <$fh>) {
	chomp($line);
	push(@{$linearray},$line);
}
close($fh);

my $rxndata;
my $rxnarray;
print "id\tnames\n";
for (my $i=0; $i < @{$linearray}; $i++) {
	my $array = [split(/\s/,$linearray->[$i])];
	print shift(@{$array})."\t";
	if (@{$array} > 0) {
		my $name = join(" ",@{$array});
		$name =~ s/;\s/|/g;
		print $name;
	}
	print "\n";
}

1;