library(optparse)

#######arguments
option_list <- list( 
    make_option(c("-i", "--input"),metavar="path", dest="otu",help="Specify the path of collapsed bacteria table",default=NULL),
    make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
    make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to compare. Required",default=NULL),
    make_option(c("-n", "--number"),metavar="int", dest="num",help="The number of most abundant species needed to be plotted, default is 20",default=15),
    make_option(c("-a", "--min-abundance"),metavar="float", dest="mina",help="The min abundance of species to be plotted, pass this will cause --number disabled",default=NULL),
    make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
    make_option(c("-b", "--by-groupMean"),metavar="logical", dest="bym",help="if T, to use the group mean to plot barplot",default=FALSE),
    make_option(c("-l", "--long-taxname"),metavar="logical", dest="long",help="if T, use the long name of species to plot heatmap",default=TRUE),
    make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files",default="./")
    )

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to relate the bacteria species and enviromental factors(numeric), and a heatmap is used to visualize the rank correlation coefficent"))

library(pheatmap)

if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}
opt$out<-paste(opt$out,"/",opt$prefix,sep="")


cluster<-ifelse(is.null(opt$map)|is.null(opt$group),TRUE,FALSE)
if(!cluster){
    map<-read.table(opt$map,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
    group<-na.omit(map[c(opt$group,"Description")])
    group<-group[order(rownames(group)),]
    group<-group[order(group[,1]),]
}


otu <- read.table(opt$otu,quote="",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")

otu<-otu[otu[,1]!="Others"&otu[,1]!="unclassified",]


#otu<-otu[!duplicated(otu[,1]),]
if(opt$long){
    rownames(otu)<-otu[,ncol(otu)]
}else{
    otu<-otu[!duplicated(otu[,1]),]
    rownames(otu)<-otu[,1]
}

#rownames(otu)<-otu[,1]

p1<-max(nchar(rownames(otu)))



otu<-t(otu[,-c(1,ncol(otu))])
sum_otu<-colSums(otu)

if(!is.null(opt$mina)){
    otu<-otu[,sum_otu>=as.numeric(opt$mina)]
}else{
    sel<-head(order(sum_otu,decreasing = T),as.numeric(opt$num))
    otu<-otu[,sel]
}

otu<-log(otu+1,base=10)

#otu<-log((otu+min(otu[otu!=0]))*10000)
#otu<-scale(otu,center=T,scale=T)

#apply(otu,2,mean)
#apply(otu,2,sd)

#print(otu)
p2<-0.5+(0.3*dim(otu)[1])+(0.06*p1)
if(cluster){
    pdf(paste(opt$out,"abundance_heatmap.pdf",sep = ""), height=ifelse(p2<50,p2,49.9),width=3+0.4*dim(otu)[2])
    pheatmap(otu,fontsize=10,border_color = "black",
             cluster_rows=T,clustering_distance_rows="euclidean",
             cluster_cols=T,clustering_distance_cols="euclidean")
    dev.off()
}else{
    otu<-otu[match(rownames(group),rownames(otu)),]
    if(opt$bym){
        otu<-apply(otu,2,function(x){tapply(x,INDEX = group[,1],mean)})
    }
    pdf(paste(opt$out,"abundance_heatmap.pdf",sep = ""), height=ifelse(p2<50,p2,49.9),width=3+0.4*dim(otu)[2])
    pheatmap(otu,fontsize=10,border_color = "black",
             cluster_cols=T,clustering_distance_cols="euclidean",cluster_rows=F)
    dev.off()
}
