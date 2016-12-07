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
for (my $i=0; $i < @{$linearray}; $i++) {
	my $array = [split(/\s/,$linearray->[$i])];
	for (my $j=0; $j < @{$array}; $j++) {
		if ($array->[$j] =~ m/^R\d\d\d$/) {
			if (defined($rxndata)) {
				push(@{$rxnarray},$rxndata);
			}
			$rxndata = {
				id => $array->[$j]
			};
		} elsif ($array->[$j] =~ m/ABAYE/) {
			$rxndata->{gpr} = $array->[$j];
			$rxndata->{gpr} =~ s/\&\&/ and /g;
			$rxndata->{gpr} =~ s/\|\|/ or /g;
		} elsif ($array->[$j] =~ m/\d+\.[\-0-9]+\.[\-0-9]+\./) {
			$rxndata->{enzyme} = [split(/\|\|\|/,$array->[$j])];
		} elsif ($array->[$j] =~ m/=/) {
			$rxndata->{equation} = $array->[$j];
		} else {
			if (defined($rxndata->{equation})) {
				push(@{$rxndata->{name}},$array->[$j]);
			} else {
				push(@{$rxndata->{pathway}},$array->[$j]);
			}
		}
	}
}

print "id\tname\tenzyme\tequation\tpathway\tgpr\n";
for (my $i=0; $i < @{$rxnarray}; $i++) {
	print $rxnarray->[$i]->{id}."\t";
	if (defined($rxnarray->[$i]->{name})) {
		print join(" ",@{$rxnarray->[$i]->{name}});
	}
	print "\t";
	if (defined($rxnarray->[$i]->{enzyme})) {
		print join("|",@{$rxnarray->[$i]->{enzyme}});
	}
	print "\t".$rxnarray->[$i]->{equation}."\t";
	if (defined($rxnarray->[$i]->{pathway})) {
		print join(" ",@{$rxnarray->[$i]->{pathway}});
	}
	print "\t";
	if (defined($rxnarray->[$i]->{gpr})) {
		print $rxnarray->[$i]->{gpr};
	}
	print "\n";
}

1;