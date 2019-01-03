#R script for lefse input file
library(optparse)
library(stringr)

option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="input",help="Specify the path of the table data file",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file",default=NULL),
    make_option(c("-g", "--category"),metavar="string",dest="group", help="Specify category name in mapping file",default="none")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "R script for generating group specific lefse input"))

filename_data<-opt$input
filename_meta<-opt$map
Grouplist<-opt$group

#setwd("~/Desktop")
#filename_data<-"all.merged.abundance.KeepID.Pathway.txt"
#filename_meta<-"ZZG.merged.mapping.txt"
#Grouplist<-'Group1'

output_temp<-gsub(pattern = "\\.txt$", "", filename_data)
#output_file="KEGG.txt"
#output_file=paste(KEGG_table_temp,"_",'Group1',,"_output.txt",sep="")

#data <- read.table(text=xx, row.name = 1, header = TRUE)
#meta <- read.table(text=yy, row.name = 1, header = TRUE , fill=TRUE, na.strings="")
data <- read.table(filename_data, row.name = 1, header = TRUE, check.names=FALSE, sep="\t")
meta <- read.table(filename_meta, row.name = 1, header = TRUE , fill=TRUE, na.strings="", check.names=FALSE, sep="\t")

dim(meta)

groups<-str_split(Grouplist,",")[[1]]
Subject <- colnames(data)
data2<-rbind(Subject,data)
rownames(data2)[1]<-'Subject'


for (group in groups){
  output_file=paste(output_temp,".",group,".lefse.txt",sep="")
  print(group)
  meta2<-meta[group]
  #print(meta3)
  meta3<-na.omit(meta2)
  #print(meta4)
  meta4<-as.data.frame(t(meta3))
  data3<-data2[,colnames(data2)%in%colnames(meta4)]
  #print(data3)
  data_final<-rbind(meta4, data3)
  #print(data_final)
  write.table(data_final,file = output_file, row.names = T, col.names = F, quote = F, sep = "\t")
}