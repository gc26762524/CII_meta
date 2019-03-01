#Sept 19, 2017. Adrian, add Mira and cdhitdup in amplification_for_task. 
#Oct 4  , 2017. Adrian, create qstat_mem_alloc_assembly
#April 9  , 2018. Adrian, modified the blast memory allocation functions to add an extra 1G when the script runs on CII cluster

our $VERSION = '1.01';

package QStatMemTracking;
use strict;
use warnings;
use POSIX;
use Math::Round qw(nearest_ceil);;
use Sys::Hostname;
use Date::Parse;
use Scalar::Util qw(looks_like_number);

my $host = hostname;

############################################################################################################################################
#
#	qstat_mem_tracking.pl : This scripts runs qstat command for a specific job id and records maximum virtual memory for every 300 seconds. 
#	
#	Test module: perl -MQStatMemTracking -e 'QStatMemTracking::qstat_mem_alloc("/home/adrian/testbench/Illumina/Stitched/PDE_HIS10_head100.extendedFrags.fastq","Cutadapt")'
#
#############################################################################################################################################

#my $human = Number::Bytes::Human->new(bs => 1024, si => 1);

#Amplification of qsub memory allocation proportional to the input file size
#(memory allocation per core, that is, the number written below will be divided by the given number of cores)
my %amplification_for_task = (
	"blastx" => "400",
	"Bowtie"  => "16",
	"Cdhit" => "6",
	"cdhitdup" => "10",
	"Cufflinks" => "16",
	"Cutadapt" => "4",
	"fastqc" => "5",
	"Filter" => "1",
	"General" => "2",
	"megablast" => "800",
	"MergeFastq" => "10",
	"Mira" => "600",
	"QC"  => "5",
	"Star" => "30",
	"Stitch" => "2",
	"Tophat2" => "4"
);
#TODO "Cdhit" => "6", is actually cdhitest, TODO replace in the name with cdhitest in all scripts that call this tool.
#TODO merge cdhit with cdhitdup (same for qstat_mem_alloc_cdhit and qstat_mem_alloc_cdhitdup). Just pay attention to modify in run_assembly.pl and generate_singleton_file
#TODO add diff value for cutadapt for SE and PE mode

#All qstat_mem_alloc are combined into the function below
sub qstat_mem_alloc {

	my $file=$_[0];
	my $task=$_[1];
	my $smp=$_[2];

	$task = "General" unless defined $task;
	$smp = 1 unless defined $smp;

	my $file_size_bytes=get_file_size_bytes($file);

	#if file size is smaller than 10MB than allocate 0.1G
	if ($file_size_bytes < 10000000){
		print ("Small input file.\n");
		print ("Memory to be allocated for this \"$task\" task: 1GB, smp: $smp\n");
		return "1G";
	}

	my $file_size_gb = ($file_size_bytes / (1024 * 1024 * 1024));
	$file_size_gb = nearest_ceil(.1, $file_size_gb);

	#if for a reason or another the file size in GB is zero than assign 0.1
	if ($file_size_gb == 0){
		$file_size_gb=0.1;
	}

	print "File size in GB: $file_size_gb (GB)\n";

	my $qstat_mem_alloc_size=$file_size_gb * $amplification_for_task{$task};

	if ( ($host =~ /tako/) && ($qstat_mem_alloc_size > 90) ){
		print "File size above the limit for tako server (Maximum amount of 90GB is allocated)\n";
		$qstat_mem_alloc_size = 90;
	}

	# mem allocation adjusted for smp job
	$qstat_mem_alloc_size=nearest_ceil(.1, $qstat_mem_alloc_size/$smp)+1;
	
	print ("Memory to be allocated for this \"$task\" task: ${qstat_mem_alloc_size}GB, smp: $smp\n");
	
	return $qstat_mem_alloc_size."G";
}

#Example call: my $MEM = QStatMemTracking::qstat_mem_alloc_assembly($in_fastq, $software, $smp, $file_size, $num_reads);
sub qstat_mem_alloc_assembly {

	my $file=$_[0];
	my $software=$_[1];
 	my $smp=$_[2];
 	my $num_reads=$_[3];

 	my $file_size_bytes=get_file_size_bytes($file);
 	my $mem="";

 	#if the file has less than 400M and less than 2M reads than assign a total min of 6G and no more than 4 smp
 	if($file_size_bytes <= 400000000 && $num_reads<2000000){
 		if ($smp > 4){
 			$smp=4;
 		}
 		$mem=(6/$smp);
 	} elsif($file_size_bytes <= 800000000 && $num_reads<4000000){
 		$mem=(8/$smp);
 	} elsif($file_size_bytes <= 1600000000 && $num_reads<8000000){
 		$mem=(16/$smp);
 	} elsif($file_size_bytes <= 3200000000 && $num_reads<16000000){
 		#File Size <= 3.2G and less than 16M reads
 		$mem=(32/$smp);
 	} else {
 		$mem=(48/$smp);
 	}	
 	
 	$mem = (nearest_ceil(.1, $mem) + 3);
 	print "Mem: ${mem}G, smp: $smp\n";
 	return ("${mem}G", $smp);
}

sub qstat_mem_alloc_bowtie {
	
	my $file=$_[0];
	my $task=$_[1];
	my $smp=$_[2];
	
	my $file_size_bytes=get_file_size_bytes($file);

	if($file_size_bytes <= 10000000){
		#minimum default memory 2G
		return "5G";
	}else{ qstat_mem_alloc($file,"Bowtie", $smp);}
}

sub qstat_mem_alloc_cdhit {
	
	my $file=$_[0];
	my $task=$_[1];
	
	my $file_size_bytes=get_file_size_bytes($file);

	if($file_size_bytes <= 10000000){
		#minimum default memory 2G
		return "2G";
	}else{ qstat_mem_alloc($file,"Cdhit");}
}

sub qstat_mem_alloc_cdhitdup {
	
	my $file=$_[0];
	my $task=$_[1];
	
	my $file_size_bytes=get_file_size_bytes($file);

	#if file size is less than 1GB then assign a default 6G for cdhitdup 
	#(cdhitdup seems to have a high mem consumption for statup)
	if($file_size_bytes <= 1000000000){
		#minimum default memory 6G
		print ("Memory to be allocated for this \"$task\" task: 8GB\n");
		return "8G";
	}else{ qstat_mem_alloc($file,"cdhitdup");}
}

sub qstat_mem_alloc_cufflinks {

        my $file=$_[0];
        my $task=$_[1];
        my $smp=$_[2];

        my $file_size_bytes=get_file_size_bytes($file);

        if($file_size_bytes <= 10000000){
                #minimum default memory 2G
                return "2G";
        }else{ qstat_mem_alloc($file,"Cufflinks", $smp);}
}

sub qstat_mem_alloc_fastqc {

	my $file=$_[0];
	my $task=$_[1];

	my $file_size_bytes=get_file_size_bytes($file);
	#if file size is less than 500MB then assign a default 1G for fastqc
	if($file_size_bytes <= 500000000){
		#minimum default memory 2G
		print ("Memory to be allocated for this fastqc task: 2GB\n");
		return "2G";
	} else {
		qstat_mem_alloc($file,$task);
	}
}

sub qstat_mem_alloc_megablast {

    my $file=$_[0];
    my $task=$_[1];
    my $smp=$_[2];

    my $mem="";

    my $file_size_bytes=get_file_size_bytes($file);

    if($file_size_bytes <= 1000000){ # files less than 1M 
        $mem="1.8";
    }elsif($file_size_bytes <= 100000000){# files less than 100M 
        $mem="2.0";
    }elsif($file_size_bytes <= 150000000){# files less than 150M 
        $mem="2.2";
    }elsif($file_size_bytes <= 200000000){# files less than 200M 
        $mem="2.4";
    }elsif($file_size_bytes <= 600000000){# files less than 600M 
        $mem="2.6";
    }elsif($file_size_bytes <= 800000000){# files less than 800M 
        $mem="2.8";
    }elsif($file_size_bytes <= 1000000000){# files less than 1G
        $mem="3.0";
    }elsif($file_size_bytes <= 2000000000){# files less than 2G
        $mem="4.0";
    }elsif($file_size_bytes <= 6000000000){# files less than 6G
        $mem="6.0";
    }else{ qstat_mem_alloc($file,"megablast", $smp);}

    #add an extra 1G for CII cluster since resources are available at no extra cost (adrian, April 9, 2018)
    if ($host =~ /cii/){
    	$mem=($mem+1);
    }

    #print "Mem: ${mem}G\n";
    return ("${mem}G");
}

sub qstat_mem_alloc_blastx {

    my $file=$_[0];
    my $task=$_[1];
    my $smp=$_[2];

    my $mem="";

    my $file_size_bytes=get_file_size_bytes($file);
        
	if($file_size_bytes <= 100000){ # files less than 100K 
		$mem="1.6";
	}elsif($file_size_bytes <= 500000){ # files less than 500K 
		$mem="1.7";
	}elsif($file_size_bytes <= 1000000){ # files less than 1M 
		$mem="1.8";
    }elsif($file_size_bytes <= 100000000){# files less than 100M 
		$mem="2.0";
    }elsif($file_size_bytes <= 200000000){# files less than 200M 
 		$mem="2.2";
    }elsif($file_size_bytes <= 600000000){# files less than 600M 
		$mem="2.6";
    }elsif($file_size_bytes <= 800000000){# files less than 800M 
		$mem="2.8";
    }elsif($file_size_bytes <= 1000000000){# files less than 1G
 		$mem="3.0";
    }else{ qstat_mem_alloc($file,"blastx", $smp);}

    #add an extra 1G for CII cluster since resources are available at no extra cost (adrian, April 9, 2018)
    if ($host =~ /cii/){
    	$mem=($mem+1);
    }

    print "Mem: ${mem}G\n";
    return ("${mem}G");
}

sub qstat_mem_alloc_merge_fastq {

	my $file=$_[0];
	my $task=$_[1];
	qstat_mem_alloc($file,"MergeFastq");
}

sub qstat_mem_alloc_prinseq {

	my $file=$_[0];
	my $task=$_[1];
	qstat_mem_alloc($file,$task);
}

sub qstat_mem_alloc_tophat2 {

    my $file=$_[0];
    my $task=$_[1];
    my $smp=$_[2];

    my $file_size_bytes=get_file_size_bytes($file);

    if($file_size_bytes <= 10000000){
            #minimum default memory 2G
            return "2G";
    }else{ qstat_mem_alloc($file,"Tophat2", $smp);}
}

sub qstat_mem_alloc_star {

	my $file=$_[0];
	my $task=$_[1];
	my $smp=$_[2];

	my $file_size_bytes=get_file_size_bytes($file);

	if($file_size_bytes <= 10000000){
	        #minimum default memory 4G
	        return "4G";
	}else{ qstat_mem_alloc($file,"Star", $smp);}
}

sub get_file_size_bytes {

	my $file=$_[0];

	if (-f $file && -e $file){
		my $file_num_bits = -s $file;
		#my $size=$human->format($file_num_bits);
        print "File size: $file_num_bits B\n";
        return $file_num_bits;
	}else{
		print "File $file does not exist\n"; exit;
	}
}

sub qstat_mem_tracking {

	my ($job_id, $outfile, $mem_requested, $smp, $time_requested) = @_;
	
	unless (defined $mem_requested) {die "Please, specify memory requested for qsub job.exit.\n"};
	unless (defined $smp) {$smp =1;}
	unless (defined $time_requested) {die "Please, specify time requested for qsub job.exit.\n"};

	my $max_v_mem = 0;
	my $job_complete = 0;
	open(OUTFILE, ">$outfile") or die "Cannot open $outfile\n";

	# Error checking for input arguments
	unless (defined $job_id){
		print OUTFILE 'Please provide the following on the command line:',"\n";
		print '	1) SGE job id',"\n"; exit;
	}

	# Log start time
	my $start_sysdate_str = localtime;
	print OUTFILE "memory usage tracking started on $start_sysdate_str\n";

	# Unless job complete flag is set to 1, loop through to check memory usage in every 300 seconds.
	do{
		my $cmd ="qstat -j $job_id";
		my $return=`$cmd 2>&1`; #capture STDOUT into a variable

		if ( $return =~ m/^Following jobs do not exist/) { #if STDOUT contains message indicating job_id no longer available
			$job_complete=1	;			   # set job_complete flag to 1
		}else{
			my @temp_line = split(/\n/,$return);      # Otherwise, parse STDOUT by line and grep the parts containing "usage" 
			my $num = scalar(@temp_line);		

			my $max_mem = 0;
			my $curr_mem = 0;
		
			for (my $i=0; $i<$num; $i++){	
				if ($temp_line[$i] =~ m/usage/){	
					my @temp_word = split(/\=/,$temp_line[$i]);	# parse out max vmem log and adjust units to "M"	
											# compare max vmem for different nodes and grep the maximum value 	
					if( $temp_word[5] =~ m/M$/){			
						$temp_word[5] =~ s/M//g;
						$curr_mem = $temp_word[5];
						if($curr_mem > $max_mem){ $max_mem = $curr_mem;}
					}elsif ( $temp_word[5] =~ m/G$/){  
						$temp_word[5] =~ s/G//g;
						$curr_mem = $temp_word[5]*1024;
						if($curr_mem > $max_mem){ $max_mem = $curr_mem;}
					}elsif( $temp_word[5] =~ m/K$/){
						$temp_word[5] =~ s/K//g;
						$curr_mem = $temp_word[5]/1024;
						if($curr_mem > $max_mem){ $max_mem = $curr_mem;}
					}else{
						$max_mem = $temp_word[5];
					}				
					
					if (looks_like_number($max_mem)) {
						if ($max_v_mem < $max_mem) { 
							$max_v_mem = $max_mem;
						}
					
					}
					print OUTFILE "$max_mem\n";	# return final max vmem from all nodes
				}
			}
			sleep 100;
		}
	}while ($job_complete==0);

	# Log end time
	my $end_sysdate_str = localtime;
	print OUTFILE "memory usage tracking ended on $end_sysdate_str\n";


	$max_v_mem = $max_v_mem/1024; # in GB instead of MB
	my $run_time_hr =str2time("$end_sysdate_str") - str2time("$start_sysdate_str");
	$run_time_hr = $run_time_hr/(60*60); # convert sec to hours
	#print $diff,"\n";
	
	
	# print expected MFCPU information
	if ($time_requested !~ m/NA/){

		my ($hh, $mm) = (0,0);
		my @time_cols=split(/\:/,$time_requested);
		
		unless($time_cols[0] eq "") {$hh = $time_cols[0];}
		unless($time_cols[1] eq "") {$mm = $time_cols[1];}
	
		# convert time format into numbers "24:: to be 24.00 hrs"		
		my $time_requested_hr = $hh + nearest_ceil(.001,$mm/60);  
		
		print_mfcpu($mem_requested, $time_requested_hr, $smp, "Expected");
	}else{
		print_mfcpu($mem_requested, $run_time_hr, $smp, "Expected");
	}

	# print observed MFCPU information	
	print_mfcpu($max_v_mem, $run_time_hr, $smp, "Observed");

	close(OUTFILE);
}

sub print_mfcpu{

	my ($mem, $time, $smp, $type) = @_;
	
	if ($type eq "Expected"){
		print OUTFILE "Expected RunTime (using qsub) :\n\n";
	
		if ($mem =~ m/G$/){	
			$mem =~ s/G//g;
			#my $total_mem_requested = $mem*$smp;		
			my $total_mem_requested = $mem;		
			$mem = sprintf("%.3f",$total_mem_requested);
		}elsif($mem =~ m/M$/){
			$mem =~ s/M//g;
			#my $total_mem_requested = ($mem/1024)*$smp;		
			my $total_mem_requested = ($mem/1024);		
			$mem = sprintf("%.3f",$total_mem_requested);
		}
	}elsif ($type eq "Observed"){

		if ($mem =~ m/G$/){ $mem =~ s/G//g;}
		elsif($mem =~ m/M$/){ $mem =~ s/M//g; $mem = $mem/1024;}
		$mem = sprintf("%.3f",$mem);
		print OUTFILE "Observed RunTime (using log file):\n\n";
	}
	
	print OUTFILE "Num. of Cores = $smp\n";
	print OUTFILE "Total Memory $type = $mem G\n";
	$time = sprintf("%.3f",$time);
	print OUTFILE "Total RunTime $type = $time hrs\n";
	
	my $mfcpu = calculate_mfcpu($mem, $smp, $time);
	
	print OUTFILE "$type MFCPU = $mfcpu\n\n";	
}

sub run_qacct_summary{
	
	my $job_id = $_[0];
	
	my $cmd ="qacct -j $job_id";
	my $return="";
	do{
		$return=`$cmd 2>&1`; #capture STDOUT into a variable
		sleep 100;
	}while($return =~ m/error/);
	
	print $return,"\n";

	my @lines = split(/\n/,$return);
	my $num_lines = scalar(@lines);
	my ($num_cores, $max_v_mem, $run_time_hr)=(0,0,0);
	
	for (my $i=0; $i<$num_lines; $i++){
		my @cols = split(/\s+/,$lines[$i]);
	
		if($lines[$i] =~ m/slots/){
			$num_cores=$cols[1];
		}elsif ($lines[$i]=~m/exit_status/){
			if ($cols[1] eq "0"){next;}
			else{print "qsub job $job_id failed while execution. exit.\n"; exit;}
		}elsif( $lines[$i] =~ m/maxvmem/){
			$max_v_mem=$cols[1];	
		}elsif($lines[$i] =~ m/ru_wallclock/){
			#$run_time_hr=nearest_ceil(.001, $cols[1]/(60*60));	
			$run_time_hr=$cols[1]/(60*60);	
		}
	}		

	return ($num_cores, $max_v_mem, $run_time_hr);
}

sub calculate_mfcpu{
	
	my ($mem,$smp, $run_time_hr)=@_;

	$mem=~s/G|M//g;
	my $mfcpu=sprintf("%.3f",($mem*$smp/3)*$run_time_hr);

	return $mfcpu;	
}

1;
