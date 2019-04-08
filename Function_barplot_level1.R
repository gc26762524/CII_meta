library(reshape2)
library(ggplot2)
library(RColorBrewer)


setwd("~/Desktop")

input <-"orthology	SS1a	SS2b	SS4d	SS6f
Cellular Processes	1043.46695971666	1130.12369575165	1127.80805731065	1121.18281700054
Genetic Information Processing	2437.41424028416	2411.48311365592	2445.06164575861	2381.29336405548
Human Diseases	1072.34340818909	1075.03490222703	1085.23570402328	1063.90694939557
Metabolism	16777.814941205	16972.4048662481	16887.013064746	17014.8696756497
Organismal Systems	745.907810475781	732.473855575721	737.041596809044	741.520367306097"

data<-read.table(text=input, header=T, row.name = 1, sep="\t")
class(data)
data2<-t(data[,order(colnames(data))])
head(data2)
class(data2)
data3<-prop.table(data2, margin=1)

rowSums(data3)
head(data3)

data4<-data3*100

data5<-melt(data4, id.var="SampleID")
colnames(data5)<-c('SampleID', 'Function', 'Proportion')
pallet<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2")[-c(7,8)],brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral"))


p<-ggplot(data5, aes(x = SampleID, y = Proportion, fill = Function)) + 
	geom_bar(stat = "identity",width = 0.7)+
  	guides(fill=guide_legend(title = NULL))+
  	scale_fill_manual(values = pallet)+
 	xlab("")+ylab("Sequence Number Percent (%)")+theme_bw()+
  	theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(),panel.border =  element_blank(),
        axis.text.x = element_text(angle = 45,size = 10,hjust = 1))+
  	scale_y_continuous(limits=c(0,101), expand = c(0, 0)
)


#ggsave(plot = p,filename = "all.merged.abundance.KeepID.Pathway.Level1.png" ,width = 10,height = 8,dpi = 300)
ggsave(plot = p,filename = "All.Function.abundance.KeepID.Pathway.Level1.pdf" ,width = 10,height = 8,dpi = 300)


