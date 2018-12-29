#/bin/sh -S

#########
#Please address any bugs to Cheng. 
#Date 2017.11.31
#########

otu_table=$1
meta_file=$2
category=$3
output_prefix=$4
taxa_filtered=$5

frequency=1000
depth=1000

if [ -z "$3" ]; then
	echo "\n\n"

	echo "Please provide following input parameters
		1) path of the OTU table file. 
		2) path of the meta data file.
		3) Group of interest from the meta file to be shown in the heatmaps.
		4) Output directory prefix name. (Required when filtered taxonomy will be specified)
		5) Taxonomy to be filtered, in the format of "tax1,tax2,tax3". (Optional)

		Sample Usage:
		sh $0 /share/data/IlluminaData/LYM/OTU_tables/Report/LYM_RTS.megablast.OTU_taxonomySummaryCounts.minHitCount_5.Bacteria.txt /share/data/IlluminaData/LYM/OTU_tables/Report/LYM_RTS_metadata.txt Group FilteredwithP_F_A Proteobacteria,Firmicutes,Actinobacteria
		"
	exit 0
else
	echo "################
	Running: sh $0 $1 $2 $3 $4 $5"
fi

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if [ -z ${4+x} ]; then 
        echo "output directory prefix is unset yet but will be auto-valued"
        output_prefix="NoPrefix"
        else
        echo "output directory prefix is set to '$4'"
fi

if [ -z ${5+x} ]; then 
	echo "taxa_filtered is unset yet but will be auto-valued"
	taxa_filtered="UNREALTAX"
	else 
	echo "taxa_filtered is set to '$5'"
fi

check_file() {
	echo "Checking file for $1 ..."
	file_name=$1
	if [ -f $file_name ]; then
		echo "File $file_name exists."
	else
		echo "File $file_name does not exist."
		exit
	fi
}

check_dir() {
	echo "Cheking the directory for $1 ..."
	dir_name=$1
	if [ ! -d "$dir_name" ]; then
		echo "Creating directory $dir_name"
		mkdir -pv $dir_name
	else
		echo "Directory $dir_name exists."
	fi
}

MAIN() {

	echo -e "\n#Enter Qiime2 enviroment"
	source activate qiime2-2018.8

	echo -e "\n#Checking files"
	check_file $otu_table
	check_file $meta_file	

	echo -e "\n#Setting up the directory structure"
	echo -e "\n#The output directory prefix is $output_prefix"
	sample_name=`echo $(basename $otu_table) | sed -e 's/\.txt$//'`
	output_dir=$(dirname $otu_table)/${sample_name}_Qiime2_output_${output_prefix}/
	check_dir $output_dir
	cp $otu_table $output_dir

	echo "##############################################################\n#Generate the figure for the percentage of annotated level"
	perl ${SCRIPTPATH}/stat_otu_tab.pl -unif min $otu_table -prefix ${output_dir}/Relative/otu_table --even ${output_dir}/Relative/otu_table.even.txt -spestat ${output_dir}/Relative/classified_stat_relative.xls
	perl ${SCRIPTPATH}/bar_diagram.pl -table ${output_dir}/Relative/classified_stat_relative.xls -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > ${output_dir}/Relative/Classified_stat_relative.svg

	echo -e "\n#Convert OTU table to biom format"
	biom convert -i $otu_table -o ${output_dir}/${sample_name}.taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

	echo -e "\n#Generate Qiime2 artifacts"
	qiime tools import   --input-path ${output_dir}/${sample_name}.taxonomy.biom --type 'FeatureTable[Frequency]'   --input-format BIOMV210Format   --output-path ${output_dir}/${sample_name}.count.qza
	qiime tools import   --input-path ${output_dir}/${sample_name}.taxonomy.biom --type "FeatureData[Taxonomy]"   --input-format BIOMV210Format   --output-path ${output_dir}/${sample_name}.taxonomy.qza

	echo -e "\n#Filter OTU table by taxonomy"
	qiime taxa filter-table   --i-table ${output_dir}/${sample_name}.count.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza --p-exclude $taxa_filtered --o-filtered-table ${output_dir}/${sample_name}.count.filtered.qza
	
	echo -e "\n#Generate barplot"
	qiime taxa barplot --i-table ${output_dir}/${sample_name}.count.filtered.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --m-metadata-file $meta_file --o-visualization ${output_dir}/${sample_name}.taxa-bar-plots.qzv
	qiime feature-table filter-features --i-table ${output_dir}/${sample_name}.count.filtered.qza   --p-min-frequency $frequency  --o-filtered-table ${output_dir}/${sample_name}.count.filtered.${frequency}.qza
	qiime taxa barplot --i-table ${output_dir}/${sample_name}.count.filtered.${frequency}.qza --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --m-metadata-file $meta_file --o-visualization ${output_dir}/${sample_name}.taxa-bar-plots.${frequency}.qzv
	
	echo -e "Conduct non-phylogenetic diversity analysis"
	qiime diversity core-metrics --i-table ${output_dir}/All.merged.OTU.count.qza  --p-sampling-depth $frequency --m-metadata-file $meta_file  --output-dir ${output_dir}/core-metrics

	echo -e "\n#Generate the heatmaps with the OTU (>= $frequency read) at different levels after collapsing."
	#for n in 2 3 4 5 6 7; do echo $n; qiime taxa collapse   --i-table ${output_dir}/${sample_name}.count.filtered.qza  --i-taxonomy ${output_dir}/${sample_name}.taxonomy.qza  --p-level $n  --o-collapsed-table ${output_dir}/${sample_name}-l${n}.qza; qiime feature-table filter-features   --i-table ${output_dir}/${sample_name}-l${n}.qza   --p-min-frequency $frequency  --o-filtered-table ${output_dir}/${sample_name}-l${n}.${frequency}.qza;  qiime feature-table heatmap --i-table ${output_dir}/${sample_name}-l${n}.${frequency}.qza --m-metadata-file $meta_file --m-metadata-column $category --o-visualization ${output_dir}/${sample_name}-l${n}.${frequency}.qzv; done;

	echo -e "\n#Generate the distance matrix and visulization artifact (rarefraction depth = $depth)."
	#qiime diversity core-metrics --i-table ${output_dir}/${sample_name}.count.filtered.qza --p-sampling-depth $depth --m-metadata-file $meta_file --output-dir ${output_dir}/diversity

	echo -e "\n#Exit Qiime2 enviroment"
	source deactivate

}

MAIN;
