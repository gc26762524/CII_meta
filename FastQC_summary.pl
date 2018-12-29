use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
Log::Log4perl->easy_init($INFO);
my $logger = Log::Log4perl::get_logger();

use lib './';
use FileUtils;

#Read more about GetOpt: http://perldoc.perl.org/Getopt/Long.html
#Testing: http://www.perl.com/pub/2005/07/14/bestpractices.html


=head1 NAME

    Generate FASTQC Summary (v0.0)

=head1 SYNOPSIS

perl FastQC_summary.pl  -f file_list

     Mandatory Options:
        -f/--fastqc-file-list file_list (absolute path for fastqc_data.txt files usually located in QC_report/SAMPLE_NAME_fastqc/fastqc_data.txt)
     Optional:

        -help            brief help message
        -man             full documentation

    Example:
    perl  perl tophat-cufflinks_warpper.pl -f file_list  -s Macaca_mulatta  
    wehre the file_list  should contain a list of "fastq files", and specise should 
    have the full pathes of bowtie 2 index files as well as a gtf file if alignment against a transcriptome is requested. 



MANDATORY parameters:

=over

=item -f/--fastq-file-list FILE

        Fastqc data list containing absolute path to fastqc_data.txt files (usually located in QC_report/SAMPLE_NAME_fastqc/fastqc_data.txt)
        NOTE: Unzip PDV_S1_L001_R2_001_fastqc.zip


=item OPTIONAL Parameters:


=item -help

        Brief help message

=item -man

        Full documentation

=back

=head1 DESCRIPTION

        Full Description
    This program will read the given xlsx file list cantaining probe names and tehir corresponding gene names


=cut

##############################
#       Declare Global Variables #
##############################
our $fastq_list_file;
our $specise;
our $fastqc;
our $RSeQC;
our $output;
our $output_nucleotide_content;
our $reference;
$output_nucleotide_content="\n";
our $delim="\t";
########################################
#       Get and Check Command Line Options #
########################################
get_options();


#########################
#   START MAIN METHOD   #
#########################
MAIN: {
  Fastqc_parser($fastq_list_file);
}#END MAIN
##############################


###############################
###       Subroutines       ###
###############################
#Get and Check Command Line Options
sub get_options{

        my $man = 0;
        my $help = 0;

        ## Parse options and print usage if there is a syntax error,
        ## or if usage was explicitly requested.
        GetOptions('help|?' => \$help, man => \$man,
                'fastq-files-list|f=s' => \$fastq_list_file);

        pod2usage(1) if $help;
        pod2usage(-verbose => 2) if $man;
        pod2usage("\nERROR: Please provide fastqc data list file \n") if not $fastq_list_file;

        if (FileUtils::file_exists($fastq_list_file) eq "no") {$logger->logdie("ERROR: file $fastq_list_file not Found!\n");}


        #print "Input files:\n";
        #print "\t fastq list        : $fastq_list_file\n";

}


sub make_PE_directory{
  my $list_file=shift;
  my @directory;
  my $sample_name;
  my @array_samples;
  my @array_sample_directories;
  if (FileUtils::file_exists($list_file) eq "no") {$logger->logdie("ERROR: file $list_file not Found!  \n");}
  open(my $fh, '<', $list_file);
  	my $file; 
          my $directory;
  	while (my $path = <$fh>) {
  	  chomp $path;
  	 ($file = $path) =~s/.*\///;
           ($directory=$path) =~m/(.*)\//;
            $directory=$1."/";
            print "directory  $1 \n";
            if($file=~m/(.+)\.fastq/)
            {
            $sample_name=$1;
            print $sample_name."--sample \n";
            push @array_samples ,$sample_name;
            push @array_sample_directories, "$directory$sample_name";
            }
  	}
  close($fh);
  return (\@array_samples,\@array_sample_directories);

}

#MAIN METHOD
sub Fastqc_parser{
my $list_file=shift;
my %summary;
my %inform;
open(my $fh, '<', $list_file);
while ( my $eachfile= <$fh>) 
{
#print "each file ".$eachfile."\n";

chomp($eachfile);
my $count_quality=0;
my $adaptor_content=0;
my $duplicate=0;
my $base_quality=0;
my $lowest_base_quality=100;
my $lowest_quality_base=100;
my $average_base_quality=0;
my $seq_length=0;
my $n_count=0;
my $Ncount=0;
my $diff_base=0;
my $nucleotide_diff;
my $nucleotide_count=0;
my $gc_content=0;
my $sum_gc_content=0;
my $read_number_GC10=0;
my $read_number_GC50=0;
my $percentage_read_number_GC10=0;
my $percentage_read_number_GC50=0;
my $highest_Ncount=0;
my $highest_Ncount_base=0;
my $Illumina_Universal=0;
my $Illumina_SmallRNA3=0;
my $Illumina_SmallRNA5=0;
my $Nextera_Transposase=0;
my $SOLID_SmallRNA=0;
my $highest_deduplication_level=0;
my $highest_deduplication_percentage=0;
my $highest_deduplication_percentage_total=0;
my $kmer=0;
my $largestKmer=0;
my $largestKmer_sequence;
my $read_quality=0;
my $read_quality_sum=0;
my $read_quality_Q20=0;
my $percentage_read_quality_Q20=0;
my $read_quality_Q30=0;
my $percentage_read_quality_Q30=0;
my $total_reads=0;
my $directory;
        ($directory=$eachfile) =~m/(.*)\//;
         $directory=$1."/";
	#$eachfile="/home/chenyu/Sample_PID-050/fastq/PID-050_CCGACAAC_BC7JPJANXX_L001_001.R1_fastqc/fastqc_data.txt";
	open(my $fh1, '<', $eachfile);
  while (my $line= <$fh1>) 
  {
	        	#print " lines : $line \n";
    {      	
			my @summary_array=split("\t",$line);
			if(defined $summary_array[1] )
			{
				chomp($summary_array[1]);
	      if(defined $summary_array[1] )
				{
					$summary{$summary_array[0]}=$summary_array[1]; #print all final filter result
				}
				if($summary_array[0]=~/>>/)
				{
					$count_quality=0;
                                        $adaptor_content=0;
					$duplicate=0;
					$n_count=0;
					$kmer=0;
					$diff_base=0;
					$read_quality=0;
					$gc_content=0;
				}
				if($count_quality==1)
				{
          if($summary_array[0]=~/\-/)
					{ 
					$base_quality+=$summary_array[1]*2;		#sum of all base quality	
					}
					else
					{
					$base_quality+=$summary_array[1];
					}
					if($lowest_base_quality>$summary_array[1])
					{
						$lowest_base_quality=$summary_array[1]; # lowest base quality
						$lowest_quality_base=$summary_array[0]; # base with lowest base quality
					}
				}
                                if($read_quality==1)
                                {
                                          $read_quality_sum+=$summary_array[1]; # sum of all read quality
					                           if($summary_array[0]<20)
					                           {
						                            $read_quality_Q20=$read_quality_sum;
					                           }	
				                          	if($summary_array[0]<30)
                                     {
                                        $read_quality_Q30=$read_quality_sum;
                                     }
                                }
                                if($adaptor_content==1)
                                {
                                          #print " adaptor content $summary_array[0] $summary_array[1] $summary_array[2] $summary_array[3] $summary_array[4] \n";
                                  				#	push(@Illumina_Universal,summary_array[0]);
                                  				#	push(@Illumina_SmallRNA3,summary_array[1]);
                                  				#	push(@Illumina_SmallRNA5,summary_array[2]);
                                  				#	push(@Nextera_Transposase,summary_array[3]);
                                  				#	push(@SOLID_SamllRNA,summary_array[4]);
                                       
					                              if($summary_array[1]>$Illumina_Universal)
                                        {
					                               	$Illumina_Universal=$summary_array[1];#location with the highest probability of Illumina_Universal adaptor
                                        }
                                        if($summary_array[2]>$Illumina_SmallRNA3)
                                        {
					                               	$Illumina_SmallRNA3=$summary_array[2];#location with the highest probability of Illumina_SmallRNA3
                                        }
                                        if($summary_array[3]>$Illumina_SmallRNA5)
                                        {
					                               	$Illumina_SmallRNA5=$summary_array[3];#location with the highest probability of Illumina_SmallRNA5
                                        }
                                        if($summary_array[4]>$Nextera_Transposase)
                                        {
				                              		$Nextera_Transposase=$summary_array[4];#location with the highest probability of Nextera_Transposase
                                        }
                                        if($summary_array[5]>$SOLID_SmallRNA)
                                        {
						                              chomp($summary_array[5]);
						                              $SOLID_SmallRNA=$summary_array[5];#location with the highest probability of SOLID_SmallRNA	
                                        }
                                }
				if($n_count==1)
				{
					$Ncount+=$summary_array[1];
					if($highest_Ncount<$summary_array[1])
					{
				   	 	$highest_Ncount=$summary_array[1];# the highest probability of N count 
              $highest_Ncount_base=$summary_array[0];# the location with highest probability of N count
					}
				} 	
				if($duplicate==1)
				{
					if($summary_array[1]>$highest_deduplication_percentage)
					{
					  $highest_deduplication_level=$summary_array[0];
					  $highest_deduplication_percentage=$summary_array[1];
						chomp($summary_array[2]);
					  $highest_deduplication_percentage_total=$summary_array[2];	
					}
				}
				if($kmer==1)
				{
					if($largestKmer< $summary_array[1])
					{
					  $largestKmer=$summary_array[1];	
			      $largestKmer_sequence=$summary_array[0];
					}
				} 
        if($diff_base==1)
        {
          my @temp=split("-",$summary_array[0]);
				  $nucleotide_count=$temp[0];
				  chomp($summary_array[4]);

          my $summary_sequence_length=$summary{"Sequence length"};
          if (index($summary_sequence_length, "-") != -1) {
            $summary_sequence_length= (split '-', $summary_sequence_length)[-1];
          } 
          #Why 16 , why -15?
					if($nucleotide_count<16 || $nucleotide_count>($summary_sequence_length-15))
          {
            $summary_array[1]=sprintf("%.2f",$summary_array[1]);
            $summary_array[2]=sprintf("%.2f",$summary_array[2]);
            $summary_array[3]=sprintf("%.2f",$summary_array[3]);
            $summary_array[4]=sprintf("%.2f",$summary_array[4]);
			      $nucleotide_diff.=$nucleotide_count.$delim.$summary_array[1].$delim.$summary_array[2].$delim.$summary_array[3].$delim.$summary_array[4].$delim;
          }
				}

                                if($gc_content==1)
                                {
                                          $sum_gc_content+=$summary_array[1];
					                               if($summary_array[0]<10)
					                               {
				                                    	$read_number_GC10=$sum_gc_content;
				                                 }
			                               		if($summary_array[0]<50)
					                              {
					                                   $read_number_GC50=$sum_gc_content;	
				                              	}

                                }
                                 if($summary_array[0] eq "#Base" && $summary_array[1] eq "Mean" )
                                {
                                  $count_quality=1;
                                }

                                if($summary_array[0] eq "#Quality" && $summary_array[1] eq "Count" )
                                {
                                  $read_quality=1;
                                }
                                 if($summary_array[0] eq "#Base" && $summary_array[1] eq "G" )
                                {
                                  $diff_base=1;
				                           $nucleotide_count=0;
			       	                     $nucleotide_diff="";
                                }
                                 if($summary_array[0] eq "#GC Content" )
                                {
                                  $gc_content=1;
                                }
                                 if($summary_array[0] eq "#Base" && $summary_array[1] eq "N-Count" )
                                {
                                  $n_count=1;
                                }

                                 if($summary_array[0] eq "#Duplication Level" && $summary_array[1] eq "Percentage of deduplicated" )
                                {
                                  $duplicate=1;
                                }
                                 if($summary_array[0] eq "#Position" && $summary_array[1] eq "Illumina Universal Adapter" )
                                {
                                  $adaptor_content=1;
                                }
                                 if($summary_array[0] eq "#Sequence" && $summary_array[1] eq "Count" )
                                {
                                  $kmer=1;
                                }
			}
		}
  }
        	close($fh1);




                $read_quality_Q20=$read_quality_sum-$read_quality_Q20;
                $read_quality_Q30=$read_quality_sum-$read_quality_Q30;
                $read_number_GC10=$sum_gc_content-$read_number_GC10;
                $read_number_GC50=$sum_gc_content-$read_number_GC50;
                $total_reads=$summary{"Total Sequences"};
		$output.=$directory; 
                $output.=$delim.$summary{"Filename"};

		my $tmp=$total_reads;
                $tmp=~s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
                $output.=$delim.$tmp;

                $output.=$delim.$summary{"Sequences flagged as poor quality"};
                $output.=$delim.$summary{"Sequence length"};
                $output.=$delim.$summary{"\%GC"};
                $output.=$delim.$summary{">>Per base sequence quality"};

    my $summary_sequence_length=$summary{"Sequence length"};
          if (index($summary_sequence_length, "-") != -1) {
            $summary_sequence_length= (split '-', $summary_sequence_length)[-1];
          } 
		$average_base_quality=$base_quality/$summary_sequence_length;	
		$average_base_quality=sprintf("%.2f",$average_base_quality);
                $output.=$delim.$average_base_quality;

                $lowest_base_quality=sprintf("%.2f",$lowest_base_quality);
                $output.=$delim.$lowest_base_quality;


                $output.=$delim.$lowest_quality_base;

                $output.=$delim.$summary{">>Per tile sequence quality"};
                $output.=$delim.$summary{">>Per sequence quality scores"};

		my $tmp1=$read_quality_Q20;
                $tmp1=~s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
                $output.=$delim.$tmp1;
                $percentage_read_quality_Q20=$read_quality_Q20/$summary{"Total Sequences"}*100;
                $percentage_read_quality_Q20=sprintf("%.2f",$percentage_read_quality_Q20);
                $output.=$delim.$percentage_read_quality_Q20;

                my $tmp2=$read_quality_Q30;
                $tmp2=~s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
                $output.=$delim.$tmp2;
                $percentage_read_quality_Q30=$read_quality_Q30/$summary{"Total Sequences"}*100;
                $percentage_read_quality_Q30=sprintf("%.2f",$percentage_read_quality_Q30);
                $output.=$delim.$percentage_read_quality_Q30;

                $output.=$delim.$summary{">>Per base sequence content"};
                $output.=$delim.$summary{">>Per sequence GC content"};

                my $tmp3=$read_number_GC10;
                $tmp3=~s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
                $output.=$delim.$tmp3;
                $percentage_read_number_GC10=$read_number_GC10/$summary{"Total Sequences"}*100;
                $percentage_read_number_GC10=sprintf("%.2f",$percentage_read_number_GC10);
                $output.=$delim.$percentage_read_number_GC10;

                my $tmp4=$read_number_GC50;
                $tmp4=~s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
                $output.=$delim.$tmp4;
                $percentage_read_number_GC50=$read_number_GC50/$summary{"Total Sequences"}*100;
                $percentage_read_number_GC50=sprintf("%.2f",$percentage_read_number_GC50);#
                $output.=$delim.$percentage_read_number_GC50;

                $output.=$delim.$summary{">>Per base N content"};
              #  $output.=$delim.$Ncount;
                $output.=$delim.$highest_Ncount;
                $output.=$delim.$highest_Ncount_base;

             #  $output.=$delim.$summary{">>Sequence Length Distribution"};


              # $output.=$delim.$summary{">>Sequence Duplication Levels"};
                $summary{"#Total Deduplicated Percentage"}=sprintf("%.2f",$summary{"#Total Deduplicated Percentage"});
                $output.=$delim.$summary{"#Total Deduplicated Percentage"};


                #$summary{"Percentage of deduplicated"}
                $highest_deduplication_percentage=sprintf("%.2f",$highest_deduplication_percentage);
                $output.=$delim.$highest_deduplication_percentage;

                $highest_deduplication_percentage_total=sprintf("%.2f",$highest_deduplication_percentage_total);
                $output.=$delim.$highest_deduplication_percentage_total;

                $output.=$delim.$summary{">>Overrepresented sequences"};
                $output.=$delim.$summary{">>Adapter Content"};

                $Illumina_Universal=sprintf("%.2f",$Illumina_Universal);
                $output.=$delim.$Illumina_Universal;

                $Illumina_SmallRNA3=sprintf("%.2f",$Illumina_SmallRNA3);
                $output.=$delim.$Illumina_SmallRNA3;

                $Illumina_SmallRNA5=sprintf("%.2f",$Illumina_SmallRNA5);
                $output.=$delim.$Illumina_SmallRNA5;

                $Nextera_Transposase=sprintf("%.2f",$Nextera_Transposase);
                $output.=$delim.$Nextera_Transposase;

                $SOLID_SmallRNA=sprintf("%.2f",$SOLID_SmallRNA);
                $output.=$delim.$SOLID_SmallRNA;

                $output.=$delim.$summary{">>Kmer Content"};
                $output.=$delim.$largestKmer;
                $output.=$delim.$largestKmer_sequence;
#	$output.="\n\n\n\n\n\nNucleotide diverse Table\n\n";
#	$output.=$nucleotide_diff;
                $output.="\n";

		$output_nucleotide_content.=$summary{"Filename"}.$delim.$nucleotide_diff."\n";;
		#print $summary{"Filename"}."\n"; 
}
        	close($fh);


          my $table_header="Directory".$delim."Filename".$delim."Total Sequences".$delim."Sequences Flagged As Poor Quality".$delim."Sequence Length".$delim."GC".$delim."Per Base Sequence Quality".$delim."Base Average Quality".$delim."Lowest Average Quality at a Position".$delim."Position with Lowest Average Base Quality".$delim."Per Tile Sequence Quality".$delim."Per Sequence Quality Scores".$delim."Read Number above Q20".$delim."Percentages of Read Number above Q20".$delim."Read Number above Q30".$delim."Percentage of Reads above Q30".$delim."Per Base Sequence Content".$delim."Per Sequence GC Content".$delim."Read Number above GC 10".$delim."Percentage Read above GC 10".$delim."Read Number above GC 50".$delim."Percentage Read above GC 50".$delim."Per Base N Content".$delim."Maximun N-Count Frequency".$delim."Position with Maximun N-Count Frequency".$delim."Total Deduplicated Percentage".$delim."Highest Deduplication Percentage".$delim."Highest Deduplication Total Percentage ".$delim."Overrepresented Sequences".$delim."Adapter_Content(%)".$delim."Illumina_Universal(%)".$delim."Illumina_SmallRNA3(%)".$delim."Illumina_SmallRNA5(%)".$delim."Nextera_Transposases(%)".$delim."OLID_SmallRNA(%)".$delim."Kmer Content".$delim."Highest Count Kmer".$delim."Highest Count Kmer Sequence\n";
          my $reference="Directory, The path for directory with fastq file
Filename, The fastq file name
Total Sequences, The total number of reads in the fastq file
Sequences Flagged As Poor Quality, If the reads are low quality
Sequence Length, The length of each read
GC, The average GC content (average out all reads)
Per Base Sequence Quality ( Quality average out all read )
Per Base Quality of Base Average Quality ( Quality average out all bases )
Lowest Average Quality at a Position, ( The lowest Base Quality along the reads )
Position with Lowest Average Base Quality,( The position with lowest Base Quality )
Per Tile Sequence Quality,( The Base Quality of reads )
Per Sequence Quality Scores,( Quality average out all reads )
Read Number above Q20,( The number of reads which is above Q20 )
Read Number above Q30,( The number of reads which is above Q30 )
Per Base Sequence Content,( The GC conetnt at base levels )
Per Sequence GC Content,( The GC conetent at read levels )
Read Number above GC 10,( The number of reads which is above GC 10% )
Read Number above GC 50,( The number of reads which is above GC 50% )
Per Base N Content,( The N content at base levels )

Maximun N-Count Frequency,( The frequency for the base has the highest N-count frequency )
Position with Maximun N-Count Frequency, ( The base position at which the N-Count Frequency is highest )


Total Deduplicated Percentage, ( The percentage of deuplication level among all reads  ) 
Highest Deduplication Level, ( Level of highest duplication )
Highest Deduplication Percentage %, ( The percentage of corresponding higest duplication level ) 
Highest Deduplication Total Percentage %, ( The percentage of highest duplication among the all reads ) 
Overrepresented Sequences, ( The sequence which are  overrepresented )
Adapter Content, ( The position of basea and their corresponding frequency for each adaptor )
Illumina Universal,
Illumina SmallRNA3,
Illumina SmallRNA5,
Nextera Transposase,
SOLID SmallRNA,
Kmer Content, ( Count of Kmers )
Ker Sequence, ( The sequence of higest frequency Kmer)
";
              my $output_nucleotide_content_header="File Name".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C".$delim."POS#".$delim."G".$delim."A".$delim."T".$delim."C";
             $output=$table_header.$output."\n".$reference."\n";
	open(my $fh2, '>', "./FastQCreport_summary.txt");	
	     print $fh2 $output."\n";
	close($fh2);
        open(my $fh3, '>', "./Nucleotide_Distribution_summary.txt");   
             print $fh3 $output_nucleotide_content_header;	
             print $fh3 $output_nucleotide_content."\n";
        close($fh3);


}

