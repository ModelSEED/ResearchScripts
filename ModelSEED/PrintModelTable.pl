use strict;
use warnings;
use DBI;
use DateTime;

$|=1;

my $db = DBI->connect("DBI:mysql:ModelDB:bio-app-authdb.mcs.anl.gov:3306","webappuser");
if (defined($db)) {
	my $models = $db->selectall_arrayref("SELECT * FROM ModelDB.MODEL", { Slice => {
		_id => 1,
		source => 1,
		status => 1,
		genome => 1,
		id => 1,
		owner => 1,
		name => 1,
		biomassReaction => 1,
		autoCompleteReactions => 1,
		autoCompleteMedia => 1,
		reactions => 1,
		associatedGenes => 1,
		gapFillReactions => 1,
		modificationDate => 1
	} });
	$db->disconnect;
	open(TABLE, "> /homes/chenry/ModelList.txt");
	for (my $i=0; $i < @{$models}; $i++) {
		my $datetime = DateTime->from_epoch(epoch => $models->[$i]->{modificationDate})->datetime();
		print TABLE $models->[$i]->{id}."\t".$models->[$i]->{genome}."\t".$models->[$i]->{owner}."\t".$models->[$i]->{status}."\t".$models->[$i]->{reactions}."\t".$models->[$i]->{biomassReaction}."\t".$models->[$i]->{gapFillReactions}."\t".$datetime."\n"; 
	}
	close(TABLE);
}

1;
