#Please address any suggestions/bugs to Adrian

use strict;
use warnings;
use Sys::Hostname;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use Log::Log4perl qw(:easy);

use lib '../modules/';
use FileUtils;

Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl::get_logger();

my $PERL="perl";

my $SCRIPT_BASEDIR=dirname $0;
print "Running $0:\n";
# VERSION

=head1 NAME

    Run cutadapt Wrapper (v1.1)

=head1 SYNOPSIS

ADD SYNOPSIS
   
Example: perl cutadapt_wrapper.pl -f samples.list -b adaptor_HiSeq.list

MANDATORY parameters:

=over

=item -f/--fastq-files-list [list.txt] 

	sample file list containing complete path names for fastq.gz files that need to be trimmed

=item -b/--adapters [adapters.list]

	List with sequence of adapter(s) ligated to the 5' or 3' end. 
	If the adapter is found within the read or overlapping the 3' end of the read, the behavior is the same as for the -a option. 
	If the adapter overlaps the 5' end (beginning of the read), the initial portion of the read matching the adapter is trimmed, but anything that follows is kept.

=item OPTIONAL Parameters:

=item -g/--adapters-5end [adapters_5end.list]

	List with sequence of adapter(s) ligated to the 5' end.

=item -a/--adapters-3end [adapters_3end.list]

	List with sequence of adapter(s) ligated to the s' end.
	
=item -O INT, --overlap=INT

	Minimum overlap length. 
	If the overlap between the read and the adapter is shorter than LENGTH, the read is not modified. 
	This reduces the number of bases trimmed purely due to short random adapter matches (default 10)

=item -d/--discard-untrimmed

	Throw away reads in which no adapter was found.

=item -m/--minimum-length N

	Throw away processed reads shorter than N bases (default 50).

=item -t/--disable-rc-trim

	Disable trimming of reverse complement adapter

=item -k/--k-mer-length

	k-mer length for computing the 5' end and 3'end k-mer distribution (default 10)

=item -help

	Brief help message

=item	-man

	Full documentation

=back

=head1 DESCRIPTION

 Trim (Illumina) adapters using cutadapt

=cut



my $overlap_len="10";
my $fastq_files_list="";
my $adapters_bothEnds_list="";
my $adapters_5end_list="";
my $adapters_3end_list="";
my $minimum_length=50;
my $disable_rc_trim=0;
my $discard_untrimmed=0;
my $kMerLength="10";
our $primer_trimmed_dir;
our $base_dir;

########################################
#	Get and Check Command Line Options #
########################################
get_options();


#########################
#	START MAIN METHOD	#
#########################
MAIN: {

	my $samples_counter=0;
	my @job_names_arr;
	my ($sample_name, $dir, $ext);
	my @exts = qw(.fa .fasta .fq .fastq .gz);

	open (SAMPLES, "<$fastq_files_list");

	while(<SAMPLES>) {	    
		chomp;
		my $this_fastq = $_;

		if(-e $this_fastq && -f $this_fastq) {

			$samples_counter++;
			print "\n~~~~ SAMPLE: $samples_counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

			($sample_name, $dir, $ext) = fileparse($this_fastq, @exts);

			print "Processing File: $sample_name\n";
			print "Working directory: $dir\n";
			print "File extenssion: $ext\n";
			
			#Get job name: Extract only filename and remove extension
		    my $job_name = (split(/\./, $sample_name))[0]."_wrapper";
		    push @job_names_arr, $job_name;

		    ($base_dir, $primer_trimmed_dir) = get_primer_trimmed_dir($this_fastq);
            my $resources;
            if (hostname =~ /hpc/) {
                	$resources="h_vmem=1G,time=23::";
            	} else {
                    $resources="h_vmem=1G";
           	}

			my $cmd = "echo perl $SCRIPT_BASEDIR/run_cutadapt.pl $this_fastq $adapters_bothEnds_list $adapters_5end_list $adapters_3end_list $overlap_len $minimum_length $disable_rc_trim $discard_untrimmed $kMerLength | qsub -V -N $job_name -l $resources -cwd -o $primer_trimmed_dir -j y";
			#print("$cmd\n");exit;
			FileUtils::run_cmd($cmd);
		}
		else{
			$logger->logdie("ERROR: Cannot find file $this_fastq");
		}

	}
	close(SAMPLES);

	print "\nFinished submitting all samples! \n";

	#Generate primer triming stats
	generate_stats(\@job_names_arr, $base_dir);

	
	print "####################################\n";
	print "# Finished processing ALL SAMPLES! #\n";
	print "####################################\n";

	print_input_parameters();
	print "Done!\n";
}


#################
#	Subroutines	#
#################

#Get and Check Command Line Options
sub get_options{

	my $man = 0;
	my $help = 0;

	## Parse options and print usage if there is a syntax error,
	## or if usage was explicitly requested.
	GetOptions('help|?' => \$help, man => \$man,
		'overlap|O=s' => \$overlap_len,
		'fastq-files-list|f=s' => \$fastq_files_list,
		'adapters|b=s' => \$adapters_bothEnds_list,
		'adapters-5end|g=s' => \$adapters_5end_list,
		'adapters-3end|a=s' => \$adapters_3end_list,
		'minimum-length|m=s' => \$minimum_length,
		'k-mer-length|k=s' => \$kMerLength,
		'discard-untrimmed|d' => \$discard_untrimmed,
		'disable-rc-trim|t' => \$disable_rc_trim) or pod2usage(2);

	pod2usage(1) if $help;
	pod2usage(-verbose => 2) if $man;

	pod2usage("\nERROR: Please provide the sample file list containing complete path names for fastq files that need to be trimmed\n") if not $fastq_files_list;

	if ( !($adapters_bothEnds_list || $adapters_5end_list || $adapters_3end_list) ){
		$logger->fatal("\nERROR: Please provide at least one list with adapters. 
		For trimming the same adapter at both 5' and 3' ends, use: -b adapters.list
		For trimming different adapters at each end, use: -g adapters_5end.list -a adapters_3end.list
		-g and -a may also be used independently\n");exit;
	}

	check_input_parameters();
}


###Check and print input parameters
sub check_input_parameters{
	
	print "Input Parameters:\n";
	print "\tInput file list                : $fastq_files_list ... checking file ... ";
	if (FileUtils::file_exists($fastq_files_list) eq "no") { $logger->logdie("ERROR: file $fastq_files_list not Found!\n");}
	else{ print "Ok\n";}

	print "\tAdapters file for both ends    : ";
	if($adapters_bothEnds_list){
		print "$adapters_bothEnds_list ... checking file ... ";
		if (FileUtils::file_exists($adapters_bothEnds_list) eq "no") {$logger->logdie("ERROR: file $adapters_bothEnds_list not Found!\n");}
		else{ print "Ok\n";}
	}
	else{
		$adapters_bothEnds_list="None";
		print "$adapters_bothEnds_list\n";
	}

	print "\tAdapters file for 5' end only  : ";
	if($adapters_5end_list){
		print "$adapters_5end_list ... checking file ... ";
		if (FileUtils::file_exists($adapters_5end_list) eq "no") {$logger->logdie("ERROR: file $adapters_5end_list not Found!\n");}
		else{ print "Ok\n";}
	}
	else{
		$adapters_5end_list="None";
		print "$adapters_5end_list\n";
	}

	print "\tAdapters file for 3' end only  : ";
	if($adapters_3end_list){
		print "$adapters_3end_list ... checking file ... ";
		if (FileUtils::file_exists($adapters_3end_list) eq "no") {$logger->logdie("ERROR: file $adapters_3end_list not Found!\n");}
		else{ print "Ok\n";}
	}
	else{
		$adapters_3end_list="None";
		print "$adapters_3end_list\n";
	}

	print "\tMinimum overlap length (nt)    : $overlap_len\n";
	print "\tMinimum read length (nt)       : $minimum_length\n";
	print "\tDiscard untrimmed reads        : ";
	if($discard_untrimmed){print "True\n";}else{print "False\n";}
	print "\tTrim reverse complement adapter: ";
	if($disable_rc_trim){print "False\n";}else{print "True\n";}
	print "\tk-mer length (to compute stats for 3' and 5' ends): $kMerLength\n";
}

#print input parameters
sub print_input_parameters{
	
	print "Input Parameters:\n";
	print "\tInput file list                : $fastq_files_list\n";
	print "\tAdapters file for 5' end only  : $adapters_5end_list\n";
	print "\tAdapters file for 3' end only  : $adapters_3end_list\n";
	print "\tMinimum overlap length (nt)    : $overlap_len\n";
	print "\tMinimum read length (nt)       : $minimum_length\n";
	print "\tDiscard untrimmed reads        : ";if($discard_untrimmed){print "True\n";}else{print "False\n";}
	print "\tTrim reverse complement adapter: ";if($disable_rc_trim){print "False\n";}else{print "True\n";}
	print "\tk-mer length (to compute stats for 3' and 5' ends): $kMerLength\n";
}


sub get_primer_trimmed_dir{

    my ($this_fastq) = @_;

	my ($this_file, $this_dir) = fileparse($this_fastq);
    my $base_dir;
    
    if ($this_dir =~ "Host_subtracted"){
		$base_dir = basename($this_dir);
    	$this_dir =~ s/$base_dir\///;
    	$base_dir = basename($this_dir);
    	$this_dir =~ s/$base_dir\///;  
    	$base_dir = basename($this_dir);
    	$this_dir =~ s/$base_dir\///;  
    }
    else{
        $base_dir = basename($this_dir);
    	$this_dir =~ s/$base_dir\///;
	}

    my $primer_trimmed_dir = "$this_dir"."Primer_trimmed";
    print "Primer trimmed output location: $primer_trimmed_dir\n";
    
    if (! -e $primer_trimmed_dir) { mkdir($primer_trimmed_dir); }
    FileUtils::check_dir($primer_trimmed_dir);
    return ($base_dir, $primer_trimmed_dir);

}

sub get_report_dir {

    my $this_dir=$primer_trimmed_dir."/";

    my $base_dir = basename($this_dir);
    $this_dir =~ s/$base_dir\///;

    $this_dir = "$this_dir/"."Report";

    return $this_dir;
}


sub generate_stats{

	my ($job_names_arr_ref, $base_dir) = @_;
	my $holdJobIDs=join(",", @$job_names_arr_ref);
	
	print "\nSubmit generate adapter statistics. Hold on the following jobs: $holdJobIDs\n";

	print "\nPrepare command for listing log files:";
	my $error_log_dir="$primer_trimmed_dir"."/error_log";
	if (! -e $error_log_dir) { mkdir($error_log_dir); }

	my $list_job_name="logs";
	my $log_files_list=${error_log_dir}."/logs.list";
	
	my $resources;
	#$resources="h_vmem=0.1G,time=:10:";
	$resources="h_vmem=0.1G";

	#print "Wait 10 seconds to allow file handlers to be closed\n";
	my $cmd="echo \"sleep 10\" | qsub -hold_jid $holdJobIDs -V -N wait10 -l $resources -cwd -o $error_log_dir/ -j y";
	my $first_waitID=FileUtils::run_cmd($cmd);
	$holdJobIDs=$holdJobIDs.",wait10";
	#print ("$first_waitID\n");
	#exit;
	#SECOND=$(qsub -W depend=afterany:$FIRST job2.pbs)
	#echo $SECOND
	#THIRD=$(qsub -W depend=afterany:$SECOND job3.pbs)
	#echo $THIRD
	
	$cmd = "echo \"ls ${error_log_dir}/*.trimmed_log.o* > $log_files_list\" | qsub -hold_jid $holdJobIDs -V -N $list_job_name -l $resources -cwd -o $error_log_dir/ -j y;";
	#print("$cmd\n");
	FileUtils::run_cmd($cmd);
	
	$holdJobIDs=$holdJobIDs.",".$list_job_name;	

	print "\nRun generate statistics:";

	$holdJobIDs.=",$list_job_name";
	my $stats_job_name="primer_trimming_stats";
	my $report_dir=get_report_dir();

	#check if directory exists. create one if not.
	FileUtils::check_dir($report_dir);	
	my $stats_out_file="$report_dir"."/$base_dir". "_primer_trimming_stats.txt";
	$cmd = "echo \"sleep 5; $PERL $SCRIPT_BASEDIR/stats/generate_adapter_statistics.pl $log_files_list > $stats_out_file\" | qsub -hold_jid $holdJobIDs -V -N $stats_job_name -l $resources -cwd -o $error_log_dir/ -j y";
	FileUtils::run_cmd($cmd);
}

