use strict;
use Data::Dumper;
use warnings;
use File::Basename;
use Sys::Hostname;
use Log::Log4perl qw(:easy);

use lib './';
use FileUtils;
use QStatMemTracking;

Log::Log4perl->easy_init($DEBUG);

my $logger = Log::Log4perl->get_logger('Starting Fastqc');
my $fastq_file = shift @ARGV;
my $smp=shift @ARGV;
$smp = 1 unless (defined $smp);

unless (defined $fastq_file) {
	print "1.Please provide a fastq file for QC analysis\n";
	exit;
}

my $fastqc_loc = FileUtils::get_exec("fastqc","./config.yml");
my $java_loc = FileUtils::get_exec("java","./config.yml");

$logger->info("Processing fastq file: $fastq_file\n");

my ($this_file, $this_dir) = fileparse($fastq_file);
my $outdir = dirname($this_dir) . "/QC_report/";

my @sample_names = split(/\./, $this_file);
my $job_name = $sample_names[0] . "_fastqc";
my $MEM=QStatMemTracking::qstat_mem_alloc_fastqc($fastq_file,"fastqc");
 
my $cmd = "echo $fastqc_loc $fastq_file -o $outdir -noextract -t $smp -d $outdir -j $java_loc | qsub -S /bin/sh -V -N $job_name -l h_vmem=$MEM -pe smp $smp -cwd -o $outdir -j y";
#print "Command: $cmd\n";

#my $log = "$outdir/$sample_names[0]" . ".fastqc" . ".log";
my $log = "$outdir/$this_file" . ".fastqc" . ".log";
$logger->info("Generating log file: $log\n");
	
my $job_id = FileUtils::run_cmd($cmd);
QStatMemTracking::qstat_mem_tracking($job_id, $log, $MEM, $smp, "NA");
