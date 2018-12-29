use strict;
use Data::Dumper;
use warnings;
use File::Basename;
use Log::Log4perl qw(:easy);

use lib './';
use FileUtils;
use QStatMemTracking;

Log::Log4perl->easy_init($DEBUG);

my $logger = Log::Log4perl->get_logger('Starting Fastqc');
my $file_list = shift @ARGV;
my $smp=shift @ARGV;
$smp = 1 unless (defined $smp);

unless (defined $file_list) {
	print "Please provide a file list containing complete path to fastq file for QC analysis\n";
	print "Please provide number of threads to be used for parallel processing (optional: no >= 2)\n";
	exit;
}

open(LIST, "<$file_list");

while(<LIST>) {
	chomp;
	my $this_fastq = $_;
	my ($this_file, $this_dir) = fileparse($this_fastq);
    	my $dir_to_process = basename($this_dir);
	my $outdir = dirname($this_dir) . "/QC_report/";
	if (! -e $outdir) { mkdir($outdir); }
	my @sample_names = split(/\./, $this_file);
	my $job_name = $sample_names[0] . "_". "$dir_to_process". "_fastqc_wrapper";
    	print $job_name, "\n";
	
	my $cmd = "echo perl run_fastqc.pl $this_fastq $smp | qsub -S /bin/sh -V -N $job_name -l h_vmem=2G -cwd -o $outdir -j y";
	my $job_id = FileUtils::run_cmd($cmd);
}

close (LIST);
