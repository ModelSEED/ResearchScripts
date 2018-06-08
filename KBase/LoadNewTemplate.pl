#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use fba_tools::fba_toolsImpl;
local $| = 1;

my $handler = fba_tools::fba_toolsImpl->new();
Bio::KBase::kbaseenv::create_context_from_client_config();
Bio::KBase::ObjectAPI::functions::set_handler($handler);

my $ws = Bio::KBase::kbaseenv::ws_client();

my $workspace = $ARGV[0];
my $template_name = $ARGV[1];

chdir("/disks/p3dev3/code/ModelSEEDDatabase/Scripts/Release/");

system("python Build_Model_Template.py ".$template_name." /disks/p3dev3/code/ModelSEEDDatabase/Templates/".$template_name);

my $template_data = Bio::KBase::ObjectAPI::utilities::FROMJSON(join("\n",@{Bio::KBase::ObjectAPI::utilities::LOADFILE("/disks/p3dev3/code/ModelSEEDDatabase/Templates/".$template_name."/".$template_name.".json")}));
$template_data->{biochemistry_ref} = "kbase/default";

my $tempobj = Bio::KBase::ObjectAPI::KBaseFBA::ModelTemplate->new($template_data);
$tempobj->parent($handler->util_store());
my $rxns = $tempobj->reactions();
for (my $i=0; $i < @{$rxns}; $i++) {
	$rxns->[$i]->compute_penalties({});
}

$handler->util_save_object($tempobj,$workspace."/".$template_name});