use strict;
use warnings;
use File::Basename;
use Sys::Hostname;
use Log::Log4perl qw(:easy);
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);

use lib '../modules/';
use FileUtils;


# VERSION
=head1 NAME

	Prinseq QC Report generation and Filtration wrapper Script

=head1 SYNOPSIS
   
	Example: perl prinseq_wrapper.pl --aim [QC/Filter] --single-end-fq-list se.prinseq.list

MANDATORY parameters:

=over

=item -a/--aim [QC/Filter]

	define the aim/purpose for prinseq run

=item -f/--single-end-fq-list [file.list]
	
	file list containing complete path names for SE fastq files

=item OPTIONAL Parameters:

=item -help

        Brief help message

=item -man

        Full documentation

=back

=head1 DESCRIPTION

Description:
Prinseq software generates QC reports and/or filter dataset with user defined criteria in "prinseq_filtration.config" file. 

=cut

#Declare global variables
my $man = 0;
my $help = 0;
my $purpose="";
my $file_list_se="";

# Parse options and print usage if there is a syntax error,
# or if usage was explicitly requested.
GetOptions('help|?' => \$help, man => \$man,
           'aim|a=s' => \$purpose,
	   'single-end-fq-list|f=s' => \$file_list_se) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

## Print usage if no arguments are given
pod2usage("\nERROR: Please provide the aim/purpose for mapping:         --aim [QC/Filter]\n") if not $purpose;
pod2usage("\nERROR: Please specify list containing complete path names for SE fastq files: --single-end-fq-list [raw_r1.list]\n") if not $file_list_se;

if ($purpose ne "QC" and $purpose ne "Filter"){
	print "\nERROR:  $purpose\n";
	print "\nERROR: Wrong value for purpose of prinseq run! Please use --aim QC or --aim Filter\n\n";
	pod2usage(1);
}

pod2usage("\nERROR: Too Many values! Check parameters!\n")  if (@ARGV != 0);

my $pwd = `pwd`;
chomp($pwd);

print "Checking file list!\n";
open (SAMPLES, "<$file_list_se") or die "ERROR: File $file_list_se not Found!\n";
my $logger = Log::Log4perl->get_logger('Starting prinseq wrapper for $file_list_se');


MAIN: {

    while(<SAMPLES>) {
    
        chomp;
        my $this_fastq = $_;
	my $filtration_config_file= "$pwd/prinseq_filtration.config";

        if(-e $this_fastq && -f $this_fastq) {
            my ($sample_file, $dir) = fileparse($this_fastq);
            my @directories = split(/Filtered*|Host_subtracted.*|Raw_fastq.*|Primer_trimmed*/, $dir);
            my $out_dir = dirname($dir);
            my $target_dir=$directories[0];
            $target_dir =~ s/\/$//g; 
            print "Target directory: $target_dir\n";
            
            if ($purpose eq "QC") {
                $out_dir = $target_dir. "/QC_report";
            } else {
                $out_dir = $target_dir. "/Filtered";
            }
            
            print "Standard output location: $out_dir\n";

            my $cmd = "cp $filtration_config_file $out_dir";
	    FileUtils::run_cmd($cmd, "General");

		my $resources;
                if (hostname =~ /hpc/) {
                        $resources="h_vmem=2G,time=23::";
                } else {
                        $resources="h_vmem=2G";
                }

            $cmd = "echo perl run_prinseq.pl $this_fastq $purpose | qsub -V -N $sample_file -l $resources -cwd -o $out_dir -e $out_dir -j y";  
		#print("$cmd\n");exit;
		FileUtils::run_cmd($cmd);
        }
        else{
            $logger->logdie("ERROR: Cannot find file $this_fastq");
        }

    }
    close(SAMPLES);
}
