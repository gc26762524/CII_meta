library(reshape2)
library(ggplot2)
library(RColorBrewer)


setwd("~/Desktop")

input <-"orthology	P1A	P1B	P2A	P2B	P3A	P3B	P4A	P4B	P5A	P5B
Cellular Processes	873.648202983645	1057.15768499962	708.328245388029	1062.14426936326	727.938344601996	863.420604335381	880.870426403297	795.159654535274	850.505643090905	717.515389209581
Genetic Information Processing	2256.09257556117	1992.95308417716	2609.9731855531	2364.47652935611	2202.1669461446	2103.17814796548	2406.44294085805	2350.00735324214	2568.96605522291	2082.69254400087
Human Diseases	859.469544312358	901.686789685348	896.795524611195	1008.83337112236	883.971896892282	868.401035151933	909.302011708902	891.900749835	939.876610328933	851.400900343522
Metabolism	13667.463856429	11642.3723302699	14718.4781940045	11304.8774368374	13976.7931462124	13542.9995912367	13620.2098947368	14227.9850086869	14019.5773510498	14193.9430732743
Organismal Systems	481.528667466912	327.844198082646	460.836986272504	341.47964535702	506.591742427993	483.90187811745	435.508183210265	455.254212885024	422.968380575738	514.166235648226
"

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


ggsave(plot = p,filename = "all.merged.abundance.KeepID.Pathway.Level1.png" ,width = 10,height = 8,dpi = 300)
ggsave(plot = p,filename = "all.merged.abundance.KeepID.Pathway.Level1.pdf" ,width = 10,height = 8,dpi = 300)


