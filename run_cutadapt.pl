use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Sys::Hostname;
use IO::Zlib;
use Memory::Usage;
use Log::Log4perl qw(:easy);

use lib './';
use FileUtils;
use QStatMemTracking;

Log::Log4perl->easy_init($DEBUG);
my $logger = Log::Log4perl->get_logger('Starting cutadapt');

#Please provide 2 files on the command line: 1) File list containing complete paths of fastq files. 2) a list of primers to be trimmed.
#Provide the following 4 arguments: 
#1) Fastq file to be trimmed.
#2) File with primers to be trimmed for both ends (None if not available)
#3) File with primers to be trimmed for 5'end (None if not available).
#4) File with primers to be trimmed for 3'end (None if not available).
#5) Overlap value for cut adapt
#6) Minimum read length (Throw away processed reads shorter than N bases.)
#7) Disable trimming of reverse complement adapter
#8) Throw away reads in which no adapter was found.
#9) k-mer length for computing the 5' end and 3'end k-mer distribution

my $this_fastq = shift @ARGV;
my $adapters_bothEnds_list = shift @ARGV;
my $adapters_5end_list = shift @ARGV;
my $adapters_3end_list = shift @ARGV;
my $overlap_len=shift @ARGV;
my $minimum_length=shift @ARGV;
my $disable_rc_trim=shift @ARGV;
my $discard_untrimmed=shift @ARGV;
my $kMerLength=shift @ARGV;

my ($primer_trimmed_dir, $report_dir, $sample_name);
my ($out_file, $info_file, $log_file);
my ($discard_untrimmed_status, $disable_rc_trim_status);

unless ((defined $this_fastq && $adapters_bothEnds_list && $adapters_5end_list && $adapters_3end_list ) && (scalar @ARGV < 4)) {
    print "Please provide at least the following parameters on command line:
    1) Fastq File to be trimmed
    2) A list of primers to be trimmed for both ends (None if not available)
    3) A list of primers to be trimmed for 5' end (None if not available)
    4) A list of primers to be trimmed for 3' end (None if not available)
    5) Overlap value for cut adapt
    6) Minimum read length (Throw away processed reads shorter than N bases.)
    7) Disable trimming of reverse complement adapter
    8) Throw away reads in which no adapter was found.
    9) k-mer length for computing the 5' end and 3'end k-mer distribution
#
    Example: perl run_cutadapt.pl sample.fastq adaptor_Illumina.list None None 10 0 0 0 10\n";
    exit;
}

$logger->info("Running $0:\n");

#my $cutadapt_loc = "/ifs/scratch/msph/cii/bl2493/Software/cutadapt-1.8.3/bin/cutadapt";
my $cutadapt_loc = FileUtils::get_exec("cutadapt","./config.yml");
print "Cutadapt_loc: $cutadapt_loc\n";

my %primers;
my %primers_order_hash;

#Check and initilaize variables
init_vars();

print("Input parameters:
\tCutadapt Path: $cutadapt_loc
\tAdapters file for both ends: $adapters_bothEnds_list
\tAdapters file for 5' end: $adapters_5end_list
\tAdapters file for 3' end: $adapters_3end_list
\tMinimum overlap length: $overlap_len
\tMinimum read length: $minimum_length
\tDiscard untrimmed reads: $discard_untrimmed_status
\tTrim reverse complement adapter: $disable_rc_trim_status
\tk-mer length (to compute stats for 3' and 5' ends): $kMerLength");

my $start_run = time();

$logger->info("Processing fastq file: $this_fastq\n");
my $mu = Memory::Usage->new();
$mu->record('Processing fastq file');

my @exts = qw(.fa .fasta .fq .fastq .gz);
my ($this_file, $this_dir, $ext) = fileparse($this_fastq, @exts);

$primer_trimmed_dir=get_primer_trimmed_dir($this_fastq);
#print "Primer trimmed output location: $primer_trimmed_dir\n"; 
if (! -e $primer_trimmed_dir) { mkdir($primer_trimmed_dir); }

my $error_log_dir="$primer_trimmed_dir"."/error_log";
#print "Primer trimmed error log output location: $error_log_dir\n";
if (! -e $error_log_dir) { mkdir($error_log_dir); }
FileUtils::check_dir($error_log_dir);

my @names = split(/\./, $this_file);
$sample_name = $names[0];
print "Sample Name:$sample_name\n";

my $lineCount = 0;

my %kMer5EndHash;
my %kMer3EndHash;
my %adaptersCountHash;

#Initialize adaptersCountHash 
foreach my $adapt_type (keys %primers_order_hash) {
    foreach my $this_primer_seq ( @{$primers_order_hash{$adapt_type}}) {
        $adaptersCountHash{$this_primer_seq}=0;
    }
}

if ($ext eq ".gz" ){
	print "Opening zip file\n";
	tie *FILE, 'IO::Zlib', "$this_fastq", "rb";
} else {
	print "Opening $ext file\n";
	open( FILE, "$this_fastq" ) || die "Could not open \n";
}

while(<FILE>){

    $lineCount++;#line 1

    my $current_seq = "";
    #get sequence
    my $line = <FILE>;$lineCount++;#line 2
    $current_seq = $line;
    $line = <FILE>;$lineCount++;#line 3
    $line = <FILE>;$lineCount++;#line 4
	#print "Current_seq: $current_seq\n";exit();
    my $fiveEndSeq=substr $current_seq, 0, $kMerLength;
    my $threeEndSeq=substr $current_seq, -$kMerLength;
    chomp $threeEndSeq;

    if ($kMer5EndHash{$fiveEndSeq}) {
        $kMer5EndHash{$fiveEndSeq}++;
    } else {
        $kMer5EndHash{$fiveEndSeq}=1;
    }

    if ($kMer3EndHash{$threeEndSeq}) {
        $kMer3EndHash{$threeEndSeq}++;
    } else {
        $kMer3EndHash{$threeEndSeq}=1;
    }

    foreach my $adapt_type (keys %primers_order_hash) {
        foreach my $this_primer_seq ( @{$primers_order_hash{$adapt_type}}) {

            my $this_primer_name = $primers{$adapt_type}{$this_primer_seq};

            if ( $current_seq=~/$this_primer_seq/ ) {
                $adaptersCountHash{$this_primer_seq}++;
            }
        }
    }
    
    #Stop processing the fastq file after 4M lines (i.e. 1M reads)
    if ( $lineCount == 4000000){
        last;
    }
}

close( FILE );


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$logger->info ("Finished parsing the fastq file.\n");
my $end_run = time();
my $run_time = $end_run - $start_run;
$logger->info("Job took $run_time seconds\n");
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$mu->record('Finished parsing the fastq file.');
$mu->dump();
print "----------------------------------------------------\n";

$logger->info ("Writing counts to files\n");
$start_run = time();

my $adaptersCountFile = "$error_log_dir/"."$sample_name".".cutadapt_counts.txt";

open(my $fh_adapt, '>', $adaptersCountFile) or die "Could not open file '$adaptersCountFile'";

print $fh_adapt "AdapterSequnce\tAdapterCount\n";
foreach my $adapterSeq (sort { $adaptersCountHash{$b} <=> $adaptersCountHash{$a} or $b cmp $a } keys %adaptersCountHash) {
    printf $fh_adapt "%-8s %s\n", $adapterSeq, $adaptersCountHash{$adapterSeq};
}

my $topCount=0;#TODO: for future choose a threshold for count number instead of printing top 30
print $fh_adapt "\n5end_kMerSequnce\tkMerCount(Top30)\n";
foreach my $kMerSeq (sort { $kMer5EndHash{$b} <=> $kMer5EndHash{$a} or $b cmp $a } keys %kMer5EndHash) {
    printf $fh_adapt "%s\t%s\n", $kMerSeq, $kMer5EndHash{$kMerSeq};
    $topCount++; if($topCount >= 30){last;}
}

$topCount=0;
print $fh_adapt "\n3end_kMerSequnce\tkMerCount(Top30)\n";
foreach my $kMerSeq (sort { $kMer3EndHash{$b} <=> $kMer3EndHash{$a} or $b cmp $a } keys %kMer3EndHash) {
    printf $fh_adapt "%s\t%s\n", $kMerSeq, $kMer3EndHash{$kMerSeq};
    $topCount++; if($topCount >= 30){last;}
}

close $fh_adapt;


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$logger->info("Finished printing stats to file: $adaptersCountFile\n");
$end_run = time();
$run_time = $end_run - $start_run;
$logger->info("Job took $run_time seconds\n");
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $counter=0;
foreach my $adapt_type (keys %primers_order_hash) {
    print "\n################################\n";
    print "Processing adapters for $adapt_type\n";

    foreach my $this_primer_seq (@{$primers_order_hash{$adapt_type}}) {

        my $this_primer_name = $primers{$adapt_type}{$this_primer_seq};
        print "\nProcessing primer $this_primer_name : $this_primer_seq\n";
        
        my $adapterCount=$adaptersCountHash{$this_primer_seq};
        print "Total count for adapter $this_primer_name: $adapterCount";

	#TODO: Convert countThresh into percentage
        my $countThresh="1000";
        if ($adapterCount >= $countThresh){

            print "\nThe adapter count is above the threshold ($countThresh)! Continue with trimming\n";
		
		#The output will be gz regardless of the input
            $out_file = "$primer_trimmed_dir/". "$sample_name".".${this_primer_name}.trimmed.fastq.gz";
            $info_file = "$primer_trimmed_dir/". "$sample_name".".${this_primer_name}.trimmed.info";
            $log_file = "$error_log_dir/". "$sample_name".".${this_primer_name}.trimmed.log";
            print "Outfile: $out_file \nInfo file: $info_file\nLog File:$log_file\n";

            my $params = "-O $overlap_len";

            if ($adapt_type eq "bothEnds"){ $params .= " -b "; }
            elsif($adapt_type eq "5end")  { $params .= " -g "; }
            elsif($adapt_type eq "5end")  { $params .= " -a "; }
            else{ die "\nERROR: Wrong adapter type\n"; }

            $params .= "${this_primer_seq}";

            if($discard_untrimmed){
                $params .= " --discard-untrimmed";
            }

            if($minimum_length){
                $params .= " --minimum-length $minimum_length";
            }

            print $params, "\n";

            my $mem=QStatMemTracking::qstat_mem_alloc($this_fastq,"Cutadapt");
            my $time="06::";
            my $job_name = "$sample_name".".${this_primer_name}".".trimmed_log";
            my $resources;
	    $resources="mem=$mem,time=$time";

            my $cmd = "echo $cutadapt_loc $params --info-file=$info_file $this_fastq -o $out_file | qsub -S /bin/sh -V -N $job_name -l $resources -cwd -o $error_log_dir -j y -sync y";

            FileUtils::run_cmd($cmd);

            print("Removing Files from Previous Steps:\n");
            system("rm -v $info_file");
            if($counter!=0){   
                system("rm -v $this_fastq");
            }
            $counter=$counter+1;
            $this_fastq = $out_file;
        }
        else{
             print "\nThe adapter count is below the threshold ($countThresh)! Continue to next adapter\n";
        }

        print "Fastq for next round of primer trimming:$this_fastq\n";

    }

}

print "\nRename the last set of files";
my $out_file_final = "$primer_trimmed_dir/${sample_name}.fq";
if ($ext eq ".gz" ){
	$out_file_final = "${out_file_final}.gz";
}

my $cmd;
if($out_file){
    $cmd = "mv $out_file $out_file_final";
}
else{
    $cmd = "cp $this_fastq $out_file_final";   
}
FileUtils::run_cmd($cmd);

print "\nDone\n";


#################
#   Subroutines #
#################

#Check and initialize variables
sub init_vars {
    $overlap_len = "10" unless defined $overlap_len;
    if ($overlap_len eq "") { $overlap_len="10";}# this is for current script only (in case it runs without wrapper

    $minimum_length = "0" unless defined $minimum_length;
    if ($minimum_length eq "") { $minimum_length="0";}

    $disable_rc_trim = 0 unless defined $disable_rc_trim;
    if ($disable_rc_trim eq "") { $disable_rc_trim=0;}
    $discard_untrimmed = 0 unless defined $discard_untrimmed;
    if ($discard_untrimmed eq "") { $discard_untrimmed=0;}

    $kMerLength = 10 unless defined $kMerLength;
    if ($kMerLength eq "" || $kMerLength eq "0") { $kMerLength=10;}

    if($discard_untrimmed){$discard_untrimmed_status="True";}else{$discard_untrimmed_status="False";}
    if($disable_rc_trim){$disable_rc_trim_status="False";}else{$disable_rc_trim_status="True";}

    if ($adapters_bothEnds_list ne "None"){

        print "Loading adapters for both ends\n";
        if (FileUtils::file_exists($adapters_bothEnds_list) eq "no") {$logger->logdie("ERROR: file $adapters_bothEnds_list not Found!\n");}
        my ($primers_bothEnds, $primer_order) = get_primer_list($adapters_bothEnds_list);
        $primers{'bothEnds'} = $primers_bothEnds;
        $primers_order_hash{'bothEnds'} = $primer_order;
    }

    if ($adapters_5end_list ne "None"){

        print "Loading 5' end adapters\n";
        if (FileUtils::file_exists($adapters_5end_list) eq "no") {$logger->logdie("ERROR: file $adapters_5end_list not Found!\n");}
        my ($primers_5End, $primer_order) = get_primer_list($adapters_5end_list);
        $primers{'5End'}=$primers_5End;
        $primers_order_hash{'5End'} = $primer_order;
    }

    if ($adapters_3end_list ne "None"){

        print "Loading 3' end adapters\n";
        if (FileUtils::file_exists($adapters_3end_list) eq "no") {$logger->logdie("ERROR: file $adapters_3end_list not Found!\n");}
        my ($primers_3End, $primer_order) = get_primer_list($adapters_3end_list);
        $primers{'3End'}=$primers_3End;
        $primers_order_hash{'3End'} = $primer_order;
    }

    if (($adapters_bothEnds_list ne "None") && ($adapters_5end_list ne "None") && ($adapters_3end_list ne "None")) {

        die "\nERROR: Please provide at least one list with adapters. 
        For trimming the same adapter at both 5' and 3' ends, use: -b adapters.list
        For trimming different adapters at each end, use: -g adapters_5end.list -a adapters_3end.list
        -g and -a may also be used independently\n";
    }
}

##################################################################################################
=head
get_primer_list 
Saves the primer sequences into an array (together with their corresponding reverse complements)
=cut
sub get_primer_list {
    
	my ($primer_list) = @_;
	open(PRIMER_LIST, "<$primer_list") || die "Can't open primer list file\n";

 	my %primers_hash;
	my ($primer_name, $primer_seq);

	my @primer_order;

    while(<PRIMER_LIST>) {
        chomp;
        next if ($_ =~ /^$/);

        
        if ($_ =~ /^>/){
            $primer_name=$_;

            #remove spaces (if any)
            $primer_name =~ s/ +//g;

            #remove the first ">" char
            $primer_name =~ s/^.//s;
        }
        else{
            $primer_seq=$_;
            $primers_hash{$primer_seq} = $primer_name;
		    push @primer_order, $primer_seq;

            if($disable_rc_trim==0){
                #Reverse complement
                my $primer_seq_rc=reverse_complement($primer_seq);
                $primers_hash{$primer_seq_rc} = ${primer_name}."_rc";
		        push @primer_order, $primer_seq_rc;
            }
        }
    }

return (\%primers_hash, \@primer_order);
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);
    
    # complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub get_primer_trimmed_dir{

    #Get passed arguments
    my ($this_fastq) = @_;

    my ($this_file, $this_dir) = fileparse($this_fastq);
    
    if ($this_dir =~ "Host_subtracted"){
        my $base_dir = basename($this_dir);
        $this_dir =~ s/$base_dir\///;
        $base_dir = basename($this_dir);
        $this_dir =~ s/$base_dir\///;  
        $base_dir = basename($this_dir);
        $this_dir =~ s/$base_dir\///;  
    }
    else{
        my $base_dir = basename($this_dir);
        $this_dir =~ s/$base_dir\///;
    }

    my $primer_trimmed_dir = "$this_dir"."Primer_trimmed";
    #print "Primer trimmed output location: $primer_trimmed_dir\n";
    
    if (! -e $primer_trimmed_dir) { mkdir($primer_trimmed_dir); }
    FileUtils::check_dir($primer_trimmed_dir);

    return $primer_trimmed_dir;

}


