use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

#################################################################
#
#	This script generates directory structure for Illumina Pipeline.
#	Usage: perl create_dir_structure.pl /share/data/IlluminaData/TEST/
#
#################################################################
my $proj=shift @ARGV;
unless (defined $proj) {
    print "Please provide a location to the project directory\n";
    exit;
}

# If project directory exists
if (-d $proj){

	print "directory $proj exists!\n";
	
	# Check if individual sub directory exists.
	check_dir("Raw_fastq");	
	check_dir("QC_report");	
	check_dir("Filtered");
	check_dir("Primer_trimmed");	
	check_dir("Assembly");	
	check_dir("Blast");	
	check_dir("Report");	
	check_dir("Spreadsheets");	
	check_dir("Kraken2");
	check_dir("AMR");
	check_dir("FMAP");

	#Copy the config.yml file for any new created project directory.
	my $config_name = "config.yml";
	my $abs_path_script = abs_path($0);
	my ($script_name, $script_dir) = fileparse($abs_path_script);
	my $abs_path_config = $script_dir. $config_name;
	print $abs_path_config, "\n";
	my $proj_spreadsheets = $proj . "Spreadsheets";
	my $CMD = "cp $abs_path_config $proj_spreadsheets";
	print $CMD, "\n";
	system($CMD);

}else{
	print "directory $proj does not exist!\n";
}

=head
	subroutine to check if directory exists & create one if it does not.		
=cut
sub check_dir{
	my ($dir_name) = @_;

	if ($proj =~ m/\/$/){ }else{$proj = $proj."/";}
	if (-d $proj.$dir_name){
		print "directory $proj$dir_name already exists!\n";
	}else{
		my $ret=system("mkdir $proj$dir_name");
		if ($ret == 0 ) {print "directory $proj$dir_name successfully created!\n";}
		else{ print "error : failed to create $proj$dir_name directory!\n";}
	}
}
