#!/usr/bin/perl

use strict;
use warnings;

my $directory = $ARGV[0];

my $clusters = {};
my $superclusters = {};

open(my $fh, "<", $directory."ClusterIDs.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$clusters->{$array->[0]} = {id => $array->[1],name => $array->[0]};
}
close($fh);

open($fh, "<", $directory."SuperClusters.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$superclusters->{$array->[0]} = {name => $array->[1], id => $array->[0], category => $array->[2]};
}
close($fh);

open($fh, "<", $directory."Cluster-SuperClusters.txt");
while (my $line = <$fh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	if (!defined($clusters->{$array->[0]})) {
		print "ERROR:Could not find ".$array->[0]."\n";
		next;
	}
	if (!defined($superclusters->{$array->[1]})) {
		print "ERROR:Could not find ".$array->[1]."\n";
		next;
	}
	$clusters->{$array->[0]}->{sc} = $superclusters->{$array->[1]};
}
close($fh);

open(my $fout, ">", $directory."Output.txt");
print $fout "Gene\tRole\tCluster\tCluster ID\tSuper cluster\tSuper cluster category\tSuper cluster ID\n";
foreach my $cluster (keys(%{$clusters})) {
	system('curl "http://aramemnon.uni-koeln.de/cluster_view.ep?cid='.$clusters->{$cluster}->{id}.'" > "'.$directory.'/temp.html"');
	open(my $fh, "<", $directory."temp.html");
	my $file = "";
	while (my $line = <$fh>) {
		chomp($line);
		$line =~ s/\t//g;
		$file .= $line;
	}
	close($fh);
	my $array = [split(/\<li\>\<a\shref\=\"javascript:openSeqSub\(\'\d+\'\);\"\stitle\=\"Sequence\"\>\<img\ssrc\=\"\.\/Gifs\/[^\.]+\.gif\"\swidth\=\"36\"\sheight\=\"20\"\salt\=\"sequence\"\/\>\<\/a\>\<a\shref\=\"javascript:openTmSub\(\'\d+\',\s0\);\"\stitle\=\"Topology\"\>\<img\ssrc\=\"\.\/Gifs\/[^\.]+\.gif\"\swidth\=\"27\"\sheight\=\"20\"\salt\=\"topology\"\/\>\&nbsp;\<\/a\>\<span\sstyle=\"font\-size:0\.9em;\"\sclass\=\"grey\"\>/,$file)];
	for(my $i=1; $i < @{$array}; $i++) {
		if ($array->[$i] =~ m/^([^\<]+)\<\/span\>\<a\sclass\=\".+\"\shref\=\"javascript:openRefGene\(\'\d+\'\);\"\stitle\=\"[^\"]+\"\>([^\<]+)\<\/a\>\<\/li\>/) {
			print $fout $1."\t".$2."\t".$clusters->{$cluster}->{name}."\t".$clusters->{$cluster}->{id}."\t".$clusters->{$cluster}->{sc}->{name}."\t".$clusters->{$cluster}->{sc}->{category}."\t".$clusters->{$cluster}->{sc}->{id}."\n";
		} elsif ($array->[$i] =~ m/^([^\<]+)\<\/span\>\<span\>([^\<]+)\<\/span\>/) {
			print $fout $1."\t".$2."\t".$clusters->{$cluster}->{name}."\t".$clusters->{$cluster}->{id}."\t".$clusters->{$cluster}->{sc}->{name}."\t".$clusters->{$cluster}->{sc}->{category}."\t".$clusters->{$cluster}->{sc}->{id}."\n";
		} else {
			print $array->[$i]."\n\n";
		}
	}
}
close($fout);