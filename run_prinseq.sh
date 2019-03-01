#/bin/sh -S

############################################################################################
# 
# Input Param:
#       $1 : File Type [fastq/fasta]
#       $2 : input file name
#       $3 : output file name
#       $4 : Full path to prinseq-lite.pl 
#       $5 : Full path to prinseq-graphs-noPCA.pl
############################################################################################

display_usage() { 

	        echo "Please provide following input parameters:
                1) input file name (R1 strand only for Paired-end sequence)
                2) execution task for prinseq run [QC/Filter]
                3) sequence data type [SE/PE]
                4) Full path to prinseq-lite.pl 
                5) Full path to prinseq-graphs-noPCA.pl

                Sample Usage:
                ./run_prinseq.sh /data/IlluminaData/PDF/Pool1/Unbiased/Primer_trimmed/GCS-076.fq.gz Filter SE /data/software/prinseq-lite-0.20.3/prinseq-lite.pl /data/software/prinseq-lite-0.20.3/prinseq-graph.pl"
        exit;
} 

if [  $# -le 4 ]; then 
	display_usage
	exit 1
fi

INFILE=$1
TASK=$2
MODE=$3
PRINSEQ_EXEC=$4
PRINSEQ_GRAPH_EXEC=$5

CONFIG_FILE=`pwd`
CONFIG_FILE=$CONFIG_FILE"/prinseq_filtration.config"
#QC_CONFIG="-graph_data -web ld,qd,gc,ns,pt,ts,aq,de,da,sc"
QC_CONFIG="-graph_data ${INFILE}.gd ld,qd,gc,ns,pt,ts,aq,de,da,sc"
SINGLETON_CONFIG="-derep 1"

PERL="perl"
FILTERED_DIR="Filtered"
gzinput="no"


run_prinseq_se(){

	echo "Running prinseq SE..."
	echo "GZ input: $gzinput; Task: $TASK"

	if [ "$gzinput" = "yes" ];then
		#gzip -dc myinputfile.fastq.gz | perl prinseq-lite.pl -verbose -fastq stdin -min_len 100 -out_good stdout -out_bad null | gzip > myoutputgood.fastq.gz
		#cmd="gzip -dc $INFILE | $PERL $PRINSEQ_EXEC -verbose -$FILETYPE stdin -min_len 100 -out_good stdout -out_bad null $CONFIG | gzip > ${OUTFILE}.gz"
		#echo $cmd
		#`$cmd`
		if [ "$TASK" = "Filter" ]; then
			gzip -dc $INFILE | $PERL $PRINSEQ_EXEC -$FILETYPE stdin -out_good stdout -out_bad null $CONFIG | gzip > ${OUTFILE}.fastq.gz
		else
			#Else TASK == QC
			echo "gzip -dc $INFILE | $PERL $PRINSEQ_EXEC -$FILETYPE stdin $OUTFILE $CONFIG"
			gzip -dc $INFILE | $PERL $PRINSEQ_EXEC -$FILETYPE stdin $OUTFILE $CONFIG	
		fi
	else
		#$PERL $PRINSEQ_EXEC -$FILETYPE $INFILE $OUTFILE $CONFIG
		cmd="$PERL $PRINSEQ_EXEC -$FILETYPE $INFILE $OUTFILE $CONFIG"
		echo $cmd
		eval $cmd
	fi
}

run_prinseq_pe(){
	$PERL $PRINSEQ_EXEC -$FILETYPE $INFILE -$FILETYPE2 $INFILE2 $OUTFILE $CONFIG
}

run_prinseq_graph(){

	if [ -e "$INFILE.gd" ]; then
		echo "Intermediate graph file $INFILE.gd exists. Proceed."
	else
		echo "Intermediate graph file $INFILE.gd does not exists. Exit."
		exit 0;
	fi

	$PERL $PRINSEQ_GRAPH_EXEC -i $INFILE.gd -html_all -o $INFILE.gd
}

read_filter_config(){

	if [ -e "$CONFIG_FILE" ]; then
		while read line; do
			FILTER_CONFIG=$FILTER_CONFIG" "$line
		done < $CONFIG_FILE
	else
		echo "Config file for Prinseq Filtration is not located on $CONFIG_FILE"
		exit 0;
	fi
	
}

derive_proj_dir(){

      echo "Input File = $INFILE"

        if [[ $INFILE == *Raw_fastq* ]]; then
                echo "Input files from Raw_fastq directory, " 
                DIR_NM=`dirname $INFILE`
                PROJ_DIR=`dirname $DIR_NM`
		elif [[ $INFILE == *Primer_trimmed* ]]; then
                echo "Input files from Primer_trimmed directory, " 
                DIR_NM=`dirname $INFILE`
                PROJ_DIR=`dirname $DIR_NM`
        elif [[ $INFILE == *Filtered* ]]; then
                echo "Input files from Filtered directory, " 
                DIR_NM=`dirname $INFILE`
                PROJ_DIR=`dirname $DIR_NM`
        elif [[ $INFILE == *Host_subtracted* ]] || [[ $INFILE == *Mapping* ]]; then
                echo "Input files from Host_subtracted/Mapping directory, "
                DIR_NM=`dirname $INFILE`
                PROJ_DIR=`dirname $(dirname $(dirname $(dirname $INFILE)))`
        elif [[ $INFILE == *Assembly* ]]; then
                echo "Input files from Assembly directory, "
                DIR_NM=`dirname $INFILE`
                PROJ_DIR=`dirname $(dirname $(dirname $INFILE))`
        else
                echo "Please, use input files from Raw_fastq/Primmer_trimmed/Host_subtracted/Mapping/Assembly directories."
                exit 0;
        fi	
	echo "Project directory:"
	echo $PROJ_DIR	
	echo "Directory containing fastq:"
	echo $DIR_NM
}

derive_file_type(){

	FILETYPE=`echo $INFILE | rev | cut -d"." -f1 | rev`

	if [ "$FILETYPE" = "gz" ]; then
		echo "gz input file"
		gzinput="yes"
                #get the second extenssion
		FILETYPE=`echo $INFILE | rev | cut -d"." -f2 | rev`
	fi
	

	if [ "$FILETYPE" = "fq" ] || [ "$FILETYPE" = "fastq" ]; then
		FILETYPE="fastq"
		echo $FILETYPE
	elif [ "$FILETYPE" = "fasta" ] || [ "$FILETYPE" = "fna" ] || [ "FILETYPE" = "fa" ]; then 
		FILETYPE="fasta"
		echo $FILETYPE
	else
		echo "Invalid file extension. Please, check input files"
		exit 0;
	fi
}

derive_sample_name(){
	echo "INFILE = $INFILE"
	SAMPLE_NM=`echo $INFILE| rev | cut -d"/" -f1 | rev`
	echo "SAMPLE_NM = $SAMPLE_NM"
	#remove extension:
	#SAMPLE_NM="${SAMPLE_NM%.*}"
	#remove second filed after dot (usualy this is the primer sequence, it may also be other name that should not be important)
	#SAMPLE_NM=`echo $SAMPLE_NM | cut -d"." -f1,3-`
	#Keep all names until the first dot (as it was initially)
	SAMPLE_NM=`echo $SAMPLE_NM | cut -d"." -f1`
	echo "SAMPLE_NM = $SAMPLE_NM"
}

derive_R2_file_name(){

        echo $INFILE

        INFILE2=`echo $INFILE|sed -e""s/R1/R2/g""`
        echo $INFILE2

        if [ -f $INFILE2 ]; then
                echo "R2 file $INFILE2 exists. Proceed in PE mode."
        else
                echo "R2 file $INFILE2 does not exist. exit."
                exit 0;
        fi
}

MAIN(){

	derive_proj_dir
	derive_file_type

	if [ "$TASK" = "QC" ] && [ "$MODE" = "SE" ]; then
		CONFIG=$QC_CONFIG
		echo $CONFIG
		OUTFILE="-out_good null -out_bad null"
	
		run_prinseq_se
		run_prinseq_graph
	
	elif [ "$TASK" = "QC" ] && [ "$MODE" = "PE" ]; then
		
		FILETYPE2=`echo $FILETYPE\2`
		CONFIG=$QC_CONFIG
		echo $CONFIG
		OUTFILE="-out_good null -out_bad null"
		
		derive_R2_file_name
		run_prinseq_pe	
		run_prinseq_graph	
	
	elif [ "$TASK" = "Filter" ] && [ "$MODE" = "SE" ]; then
		
		derive_sample_name
		read_filter_config
		
		CONFIG=$FILTER_CONFIG
		echo $CONFIG
		
		if [ "$gzinput" = "yes" ];then
        		OUTFILE="$PROJ_DIR/$FILTERED_DIR/$SAMPLE_NM.good"
		else
			OUTFILE="-out_good $PROJ_DIR/$FILTERED_DIR/$SAMPLE_NM.good -out_bad null"    
        	fi

		run_prinseq_se
	
	elif [ "$TASK" = "Filter" ] && [ "$MODE" = "PE" ] ; then
		
		derive_sample_name
		read_filter_config
		
		FILETYPE2=`echo $FILETYPE\2`
		CONFIG=$FILTER_CONFIG
		echo $CONFIG
		#for PE data host subtraction is after filtration so expects a different directory structure.
		OUTFILE="-out_good $PROJ_DIR/$FILTERED_DIR/$SAMPLE_NM.good -out_bad $PROJ_DIR/$FILTERED_DIR/$SAMPLE_NM.bad"
		echo $OUTFILE
		derive_R2_file_name
		run_prinseq_pe

        elif [ "$TASK" = "Singletons" ]; then

                derive_sample_name
                CONFIG=$SINGLETON_CONFIG
		echo $CONFIG
                OUTFILE="-out_good $DIR_NM/$SAMPLE_NM.singletons -out_bad null"
		echo $OUTFILE
                run_prinseq_se
	else
		echo "Inappropriate Task and/or mode selected for prinseq run. Please, try QC or Filter task and specify data type whether it is SE or PE data."
		exit 0;
	fi
}

MAIN;
