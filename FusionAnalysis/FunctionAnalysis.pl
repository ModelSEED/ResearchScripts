#!/usr/bin/perl

$|=1;
my $directory = $ARGV[0];
my $genomedir = $ARGV[1];

my $array;
open(my $fh, "<", $genomedir."GenomeList.txt");
while (my $line = <$fh>) {
	chomp($line);
	push(@{$array},$line);
}
close($fh);

my $functionhash;
my $FunctionGenomes;
for (my $i=0; $i < @{$array}; $i++) {
	print "Loading ".$i.":".$array->[$i]."\n";
	open(my $fhh, "<", $genomedir.$array->[$i]);
	my $line = <$fhh>;
	my $GeneCDDs = {};
	my $newcdds = {};
	while (my $line = <$fhh>) {
		chomp($line);
		my $items = [split(/\t/,$line)];
		$items->[6] =~ s/\s*#.+//;
		my $fnarray = [sort(split(/\s*;\s+|\s+[\@\/]\s+/,$items->[6]))];
		for (my $j=0; $j < @{$fnarray}; $j++) {
			$fnarray->[$j] = lc($fnarray->[$j]);
			$fnarray->[$j] =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
			$fnarray->[$j] =~ s/\s//g;
			$fnarray->[$j] =~ s/\#.*$//g;
			$fnarray->[$j] =~ s/\(ec\)//g;
			if (!defined($functionhash->{$fnarray->[$j]})) {
				$functionhash->{$fnarray->[$j]} = [0,0];
			}
			if (!defined($FunctionGenomes->{$fnarray->[$j]}->{$array->[$i]})) {
				$FunctionGenomes->{$fnarray->[$j]}->{$array->[$i]} = 1;
				$functionhash->{$fnarray->[$j]}->[1]++;
			}
			$functionhash->{$fnarray->[$j]}->[0]++;
		}
	}
}

open ( my $out, ">", $directory."/FunctionAnalysis.txt");
print $out "Function\tCount\tGenome count\n";
foreach my $function (keys(%{$functionhash})) {
	print $out $function."\t".$functionhash->{$function}->[0]."\t".$functionhash->{$function}->[1]."\n";
}
close($out);

1;