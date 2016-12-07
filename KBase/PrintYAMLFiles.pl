#!/usr/bin/perl -w

use strict;
use JSON;
use Data::Dumper;

my $dropdownoptions = {
	build_metabolic_model => {
		template_id => {
			options => [
				{
		          "value" => "auto",
		          "display" => "Automatic selection",
		          "id" => "auto",
		          "ui_name" => "Automatic selection"
		        },
		        {
		          "value" => "gramneg",
		          "display" => "Gram negative",
		          "id" => "gramneg",
		          "ui_name" => "Gram negative"
		        },
		        {
		          "value" => "grampos",
		          "display" => "Gram positive",
		          "id" => "grampos",
		          "ui_name" => "Gram positive"
		        },
		        {
		          "value" => "plant",
		          "display" => "Plant",
		          "id" => "plant",
		          "ui_name" => "Plant"
		        }
			]
		}
	},
	propagate_model_to_new_genome => {
		translation_policy => {
			options => [
				{
		          "value" => "add_reactions_for_unique_genes",
		          "display" => "Add unique gene reactions",
		          "id" => "add_reactions_for_unique_genes",
		          "ui_name" => "Add unique gene reactions"
		        },
		        {
		          "value" => "translate_only",
		          "display" => "Only translate model",
		          "id" => "translate_only",
		          "ui_name" => "Only translate model"
		        },
		        {
		          "value" => "reconcile",
		          "display" => "Merge in KBase model",
		          "id" => "reconcile",
		          "ui_name" => "Merge in KBase model"
		        }
			]
		}
	}
};

my $methodcategories = {
	run_flux_balance_analysis => ["active","metabolic_modeling"],
	propagate_model_to_new_genome => ["active","metabolic_modeling","comparative_genomics"],
	build_metabolic_model => ["active","metabolic_modeling"],
	simulate_growth_on_phenotype_data => ["active","metabolic_modeling"],
	compare_fba_solutions => ["active","metabolic_modeling"],
	gapfill_metabolic_model => ["active","metabolic_modeling"],
	merge_metabolic_models_into_community_model => ["active","metabolic_modeling","communities"],
	compare_flux_with_expression => ["active","metabolic_modeling","expression"],
	check_model_mass_balance => ["active","metabolic_modeling"],
	edit_metaboli_model => ["active","metabolic_modeling"],
	edit_or_create_media => ["active","metabolic_modeling"]
};

my $subselectionparams = {
	target_reaction => {
		"subdata_selection"=> {
	        "parameter_id" => "fbamodel_id",
			"subdata_included" => ["modelreactions/[*]/id", "modelreactions/[*]/name"],
			"path_to_subdata"=> ["modelreactions"],
			"selection_id" => "id",
			"selection_description" => ["name"],
			"description_template" =>"- {{name}}",
			"additional_options" => ["bio1 - Biomass"],
		},
		"multiselection"=>"false",
		"show_src_obj"=>"true",
		"allow_custom"=>"true"
	},
	media_supplement_list => {
		"subdata_selection"=> {
	        "parameter_id" => "fbamodel_id",
			"subdata_included" => ["modelcompounds/[*]/id", "modelcompounds/[*]/name"],
			"path_to_subdata"=> ["modelcompounds"],
			"selection_id" => "id",
			"selection_description" => ["name"],
			"description_template" =>"- {{name}}"
		},
		"multiselection"=>"false",
		"show_src_obj"=>"true",
		"allow_custom"=>"false"
	},
	expression_condition => {
		"subdata_selection"=> {
	        "parameter_id" => "expseries_id",
			"subdata_included" => ["data/col_ids"],
			"path_to_subdata"=> ["data","col_ids"],
			"selection_id" => "id"
		},
		"multiselection"=>"false",
		"show_src_obj"=>"true",
		"allow_custom"=>"false"
	},
#	feature_ko_list => 
	reaction_ko_list => {
		"subdata_selection"=> {
	        "parameter_id" => "fbamodel_id",
			"subdata_included" => ["modelreactions/[*]/id", "modelreactions/[*]/name"],
			"path_to_subdata"=> ["modelreactions"],
			"selection_id" => "id",
			"selection_description" => ["name"],
			"description_template" =>"- {{name}}"
		},
		"multiselection"=>"true",
		"show_src_obj"=>"true",
		"allow_custom"=>"false"
	}
};

#Read methods
my $parameters;
my $jsonfiles;
open(my $fhh, "<", "/Users/chenry/code/fba_tools/ui/narrative/Parameters.txt");
my $line = <$fhh>;
while ($line = <$fhh>) {
	chomp($line);
	my $array = [split(/\t/,$line)];
	$parameters->{$array->[0]}->{$array->[1]} = {
		id => $array->[1],
		method => $array->[0],
		"default" => $array->[2],
		optional => $array->[3],
		advanced => $array->[4],
		array => $array->[5],
		type => $array->[6],
		output => $array->[7],
		name => $array->[8],
		description => $array->[10],
		placeholder => $array->[9]
	};
	if (!defined($jsonfiles->{$array->[0]})) {
		$jsonfiles->{$array->[0]} = {
			name => $array->[0],
			ver => "1.0.0",
			authors => ["chenry"],
			contact => "help\@kbase.us",
			visible => "true",
			categories => $methodcategories->{$array->[0]},
			widgets => {
				input => undef,
				output => "kbaseTabTable"
			},
			job_id_output_field => "docker",
			parameters => [],
			behavior => {
				"service-mapping" => {
					url => "",
				    name => "fba_tools",
				    method => $array->[0],
				    input_mapping => [{
				          narrative_system_variable => "workspace",
				          target_property => "workspace"
				   	}],
				    output_mapping => [{
				          narrative_system_variable => "workspace",
				          target_property => "ws"
				    }],
				}
			}
		};
		if ($array->[0] eq "compare_flux_with_expression") {
			$jsonfiles->{$array->[0]}->{widgets}->{output} = "kbaseExpressionAnalysis";
		} elsif ($array->[0] eq "check_model_mass_balance") {
			$jsonfiles->{$array->[0]}->{widgets}->{output} = "kbaseReportView";
			$jsonfiles->{$array->[0]}->{behavior}->{"service-mapping"}->{output_mapping}->[0]->{target_property} = "workspace_name";
			$jsonfiles->{$array->[0]}->{behavior}->{"service-mapping"}->{output_mapping}->[1]->{service_method_output_path} = [0,"report_name"];
			$jsonfiles->{$array->[0]}->{behavior}->{"service-mapping"}->{output_mapping}->[1]->{target_property} = "report_name";
		}
	}
	my $paramobj = {
		id => $array->[1],
    	optional => "false",
    	advanced => "false",
    	allow_multiple => "false",
    	default_values => [ $parameters->{$array->[0]}->{$array->[1]}->{"default"} ],
    	field_type => "text",
	    text_options => {
	      valid_ws_types => []
	    }
	};
	if (defined($subselectionparams->{$array->[1]})) {
		$paramobj->{textsubdata_options} = $subselectionparams->{$array->[1]};
		$paramobj->{field_type} = "textsubdata";
	}
	if ($parameters->{$array->[0]}->{$array->[1]}->{optional} eq "1") {
		$paramobj->{optional} = "true";
	}
	if ($parameters->{$array->[0]}->{$array->[1]}->{advanced} == 1) {
		$paramobj->{advanced} = "true";
	}
	if ($parameters->{$array->[0]}->{$array->[1]}->{array} == 1) {
		$paramobj->{allow_multiple} = "true";	
	}
	if ($parameters->{$array->[0]}->{$array->[1]}->{type} eq "checkbox") {
		$paramobj->{field_type} = "checkbox";
		$paramobj->{checkbox_options} = {
			checked_value => 1,
			unchecked_value => 0 
		};
	}
	if ($parameters->{$array->[0]}->{$array->[1]}->{type} =~ m/dropdown/) {
		$paramobj->{field_type} = "dropdown";
		$paramobj->{dropdown_options} = $dropdownoptions->{$array->[0]}->{$array->[1]};	
	} elsif ($parameters->{$array->[0]}->{$array->[1]}->{type} =~ m/\./) {
		$paramobj->{text_options} = {
			valid_ws_types => [$parameters->{$array->[0]}->{$array->[1]}->{type}], 
		};
		if ($parameters->{$array->[0]}->{$array->[1]}->{output} == 1) {
			$paramobj->{text_options}->{is_output_name} = "true";
			push(@{$jsonfiles->{$array->[0]}->{behavior}->{'service-mapping'}->{output_mapping}},{
				constant_value => $parameters->{$array->[0]}->{$array->[1]}->{type},
		        target_property => "type"
			});
		    push(@{$jsonfiles->{$array->[0]}->{behavior}->{'service-mapping'}->{output_mapping}},{
				input_parameter => $array->[1],
		        target_property => "obj"
			});
		}
    } elsif ($parameters->{$array->[0]}->{$array->[1]}->{type} eq "float" || $parameters->{$array->[0]}->{$array->[1]}->{type} eq "int") {
    	$paramobj->{text_options} = {validate_as => $parameters->{$array->[0]}->{$array->[1]}->{type}};
    }
	push(@{$jsonfiles->{$array->[0]}->{parameters}},$paramobj);
	push(@{$jsonfiles->{$array->[0]}->{behavior}->{'service-mapping'}->{input_mapping}},{
		input_parameter => $array->[1],
        target_property => $array->[1]
	});
}
close($fhh);
open(my $fh, "<", "/Users/chenry/code/fba_tools/ui/narrative/MethodsData.txt");
my $currentmethod;
my $onprefix = 1;
my $prefix;
my $suffix;
while (my $line = <$fh>) {
	if ($line =~ m/\{(.+)\}/) {
		if (defined($currentmethod)) {
			&print_yampl($currentmethod,$prefix,$suffix,$parameters->{$currentmethod},$jsonfiles->{$currentmethod});
		}
		$currentmethod = $1;
		$prefix = "";
		$suffix = "";
		$onprefix = 1;
	} elsif ($line =~ m/---/) {
		$onprefix = 0;
	} elsif ($onprefix == 1) {
		$prefix .= $line;
	} elsif ($onprefix == 0) {
		$suffix .= $line;
	}
}

sub print_yampl {
	my ($fmethod,$fprefix,$fsuffix,$fparameters,$jsonfiles,$icon) = @_;
	print "2:".$fmethod."\n";
	my $JSON = JSON::XS->new->ascii->pretty;
	open(my $fn, ">", "/Users/chenry/code/fba_tools/ui/narrative/methods/".$fmethod."/spec.json");
	my $jsontext = $JSON->encode($jsonfiles);
	$jsontext =~ s/\"true\"/true/g;
	$jsontext =~ s/\"false\"/false/g;
	print $fn $jsontext;
	close($fn);	
	open(my $fo, ">", "/Users/chenry/code/fba_tools/ui/narrative/methods/".$fmethod."/display.yaml");
	print $fo $fprefix;
	print $fo "parameters :\n";
	foreach my $param (keys(%{$fparameters})) {
    	print $fo "    ".$param." :\n";
    	print $fo "        ui-name : |\n";
    	print $fo "            ".$fparameters->{$param}->{name}."\n";
    	print $fo "        short-hint : |\n";
    	print $fo "            ".$fparameters->{$param}->{description}."\n";
    	print $fo "        long-hint : |\n";
    	print $fo "            ".$fparameters->{$param}->{description}."\n";
		if (defined($fparameters->{$param}->{placeholder}) && length($fparameters->{$param}->{placeholder}) > 0) {
			print $fo "        placeholder : |\n";
    		print $fo "            ".$fparameters->{$param}->{placeholder}."\n";
		}
		print $fo "\n\n";
	}
    print $fo $fsuffix;
	close($fo);
}