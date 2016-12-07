#!/usr/bin/perl

use strict;
use warnings;
use Bio::KBase::workspace::ScriptHelpers qw(printObjectInfo get_ws_client workspace workspaceURL parseObjectMeta parseWorkspaceMeta printObjectMeta);
use Bio::KBase::fbaModelServices::ScriptHelpers qw(load_file get_ws_objects_list fbaws get_fba_client runFBACommand universalFBAScriptCode );

my $filename = $ARGV[0];
my $data = load_file($filename);
my $headings = [qw(
ID
Created
Updated
Name
CASS
Class
Metabolism
Organism
Target
ActionType
Sequence
)];

print join("\t",@{$headings})."\n";
my $drugdata;
for (my $i=0; $i < @{$data}; $i++) {
	my $line = $data->[$i];
	if ($line =~ m/\<drug.+created=\"(.+)\".+updated=\"(.+)\"/) {
		my $create = $1;
		my $update = $2;
		if (defined($drugdata->{ID}) && $drugdata->{Organism} =~ m/[bB]acteria/) {
			for (my $j=0; $j < @{$headings}; $j++) {
				if ($j > 0) {
					print "\t";
				}
				if (defined($drugdata->{$headings->[$j]})) {
					print $drugdata->{$headings->[$j]};
				}
			}
			print "\n";
		}
		$drugdata = {
			Created => $create,
			Updated => $update
			
		};
	} elsif ($line =~ m/\<affected-organism\>(.+)\<\/affected-organism\>/) {
		$drugdata->{Organism} = $1;
	} elsif ($line =~ m/\<drugbank-id\sprimary=\"true\"\>(DB\d+)\<\/drugbank-id\>/) {
		$drugdata->{ID} = $1;
	} elsif ($line =~ m/\<cas-number\>(.+)\<\/cas-number\>/) {
		$drugdata->{CASS} = $1;
	} elsif ($line =~ m/\<name\>(.+)\<\/name\>/) {
		$drugdata->{Name} = $1;
	} elsif ($line =~ m/\<kingdom\>(.+)\<\/kingdom\>/) {
		$drugdata->{Class} = $1.";";
	} elsif ($line =~ m/\<superclass\>(.+)\<\/superclass\>/) {
		$drugdata->{Class} .= $1.";";
	} elsif ($line =~ m/\<class\>(.+)\<\/class\>/) {
		$drugdata->{Class} .= $1.";";
	} elsif ($line =~ m/\<subclass\>(.+)\<\/subclass\>/) {
		$drugdata->{Class} .= $1;
	} elsif ($line =~ m/\<metabolism\>(.+)\<\/metabolism\>/) {
		$drugdata->{Metabolism} = $1;
	} elsif ($line =~ m/\<target[^\/]*\>/) {
		while ($line !~ m/\<\/target/ && defined($line)) {
			$i++;
			$line = $data->[$i];
			if ($line =~ m/amino-acid-sequence\sformat=\"FASTA\"\>\&gt;(.+)/) {
				$drugdata->{Target} .= $1;
				$drugdata->{Sequence} = "";
				while ($line !~ m/(.+)\<\/amino-acid-sequence\>/ && defined($line)) {
					$i++;
					$line = $data->[$i];
					if ($line =~ m/(.+)\<\/amino-acid-sequence\>/) {
						$drugdata->{Sequence} .= $1;
					} else {
						$drugdata->{Sequence} .= $line;
					}
				}
			} elsif ($line =~ m/\<name\>(.+)\<\/name\>/) {
				$drugdata->{Target} .= $1;
			} elsif ($line =~ m/\<action\>(.+)\<\/action\>/) {
				$drugdata->{ActionType} .= $1;
			}
		}
	}
}