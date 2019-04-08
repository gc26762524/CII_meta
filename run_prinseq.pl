use strict;
use warnings;
#use Bio::SeqIO;
use Data::Dumper;
use File::Basename;
use Sys::Hostname;
use Log::Log4perl qw(:easy);
use Cwd;

use lib './';
use FileUtils;
use QStatMemTracking;

my $PERL= "perl";

Log::Log4perl->easy_init($DEBUG);
###################################################################################################################
#
#       This is a master script for running prinseq either for QC report generation or prinseq Filteration
#
#	Usage: 
#	echo "/nfs/apps/perl/current/bin/perl prinseq_wrapper.pl filter.list Filter" | qsub -l mem=1G,time=01:: -cwd -o /ifs/scratch/msph/cii/bl2493/Projects/Illumina/TEST/Filtered/ -e /ifs/scratch/msph/cii/bl2493/Projects/Illumina/TEST/Filtered/ -V -N filter_test
###################################################################################################################

# Input Arguments as Global variables
my $fastq= shift @ARGV;
my $task = shift @ARGV;
my $mode = "SE";

# Error checking for input arguments
unless (defined $fastq && $task) { 
	print 'Please provide the following on the command line: 
		1) Complete path to the input fastq file
		2) Execution task [QC/Filter]',"\n";
       	exit;
}

print "Processing $fastq for $task in $mode mode\n";



MAIN();

sub MAIN{
	
	print "Fastq file = $fastq\n";  #get a fastq
	my ($file, $directory) = fileparse($fastq);  # get the location of fastq
	my $job_type  = basename($directory);
	my @directories = split(/Filtered|Host_subtracted|Raw_fastq|Primer_trimmed/, $directory);  # Remove the last dirname from the current path
    
	#print Dumper @directories;
	my @names = split(/\./, $file);    #  get the sample name by splitting the file name
	my $target_dir=$directories[0];    #   get the location of the project directory
   	$target_dir =~ s/\/$//g; 
	print "Target directory = $target_dir\n";
	print $job_type, "\n";
	my $job_name=$names[0]."_"."$job_type"."_".$task;  #   job name is the sample name and the task type
	print $job_name, "\n";
	my $mem=QStatMemTracking::qstat_mem_alloc_prinseq($fastq,$task);    # determine the memory for the job based on the file size.
	my $dir = "";
   	my $smp = 1; # prinseq does not allow smp jobs. $smp variable defaulted to 1.

 
	if ($task eq "QC") {
		$dir = "$target_dir/"."QC_report"; # determine the location of the QC reports.
	} else {
		$dir = "$target_dir/"."Filtered";  # determine the location of the Filtered the output.
	}
	print "Directory = $dir\n";
    
	# Run run_prinseq.sh script for given task (QC or Filter)

	my $resources;
        if (hostname =~ /hpc/) {
                $resources="h_vmem=$mem,time=:10:";
        } else {
                $resources="h_vmem=$mem";
        }

	print "Running Prinseq,\n";

	#my $PRINSEQ_EXEC=FileUtils::get_exec("prinseq","./config.yml");
	#my $PRINSEQ_GRAPH_EXEC=FileUtils::get_exec("prinseq_graph","./config.yml");
	my $PRINSEQ_EXEC="/home/cheng/softwares/miniconda2/bin/prinseq-lite.pl";
        my $PRINSEQ_GRAPH_EXEC="/home/cheng/softwares/miniconda2/bin/prinseq-graphs.pl";
	my $cmd = "qsub -S /bin/sh \-V \-N $job_name \-l $resources \-cwd \-o $dir \-e $dir -j y run_prinseq.sh $fastq $task $mode $PRINSEQ_EXEC $PRINSEQ_GRAPH_EXEC";
	#print "$cmd\n";exit;
	my $job_id = FileUtils::run_cmd($cmd);
    
	if ($task eq "QC"){

		# Writing log file	
		my $log = $dir. "/". $names[0].".QC.log";
		QStatMemTracking::qstat_mem_tracking($job_id, $log, $mem, $smp, "NA");
			
		$target_dir = $dir;
		FileUtils::check_dir($target_dir);			
	
		#move files to designated directories
		FileUtils::move_files("$dir$names[0]*.QC.log", $target_dir);
		FileUtils::move_files("$fastq.gd", $target_dir);
		FileUtils::move_files("$fastq.gd.html", $target_dir);
	
	}elsif($task eq "Filter"){
		
		# Writing log file	
		my $log = $dir. "/".$names[0].".Filter.log";	
		QStatMemTracking::qstat_mem_tracking($job_id, $log, $mem, $smp, "NA");
		
		$target_dir = $dir;
		FileUtils::check_dir($target_dir);			
		FileUtils::move_files("$dir/$names[0]*.good*fastq.gz", $target_dir);

		#move files to designated directories
		FileUtils::move_files("$dir/$names[0]*.Filter.log", $target_dir);
	}
	#move files to designated directories
	FileUtils::move_files("$dir/$names[0]*.o$job_id", $target_dir);
} 
