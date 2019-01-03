library(reshape2)
library(ggplot2)
library(RColorBrewer)


setwd("~/Desktop")

input <-"orthology	C10	M11	D10	C12	M13	D11	C13	M6	D2	C14	M7	D6	C15	M8	D8	C9	M9	D9
Cellular Processes	861.240037842928	900.238851397065	789.986363076506	800.776326285526	780.050520633442	722.064897238334	894.879248788463	831.831422434276	1210.8464182072	833.258412685652	857.459865127599	967.629810669541	778.650821151109	835.175016098921	1141.1490161608	931.911422443146	949.145622719107	831.107731147403
Genetic Information Processing	3056.99537905976	3213.81927518111	2998.92333699773	3182.76913138293	2961.33229979732	2993.74010660231	2975.97816452623	2801.68015685603	3112.70037952464	2910.36076213753	2880.84393001708	3033.59186301859	2851.22837762121	2977.59413810497	3063.35755216125	2900.93507985528	2673.52314278072	3024.59949091225
Human Diseases	894.46167675502	936.006975977204	927.76451272383	873.508037359363	957.865160151697	925.36912493217	861.687060806168	951.426108752615	881.328702661238	917.665968708819	948.612053568706	918.149988191629	978.111541431084	948.273287183796	869.691721296139	944.044074286991	915.71358262966	928.248660725993
Metabolism	13973.3267038654	13919.4460688275	13717.9911650167	13389.3256418622	14457.9639373768	14482.6619424078	13738.076699886	13942.4618436512	13356.8772734065	15038.1042336494	14587.0233407459	13400.4671048778	15599.7537742368	14529.9756899752	13444.1992991185	14948.8413422579	12959.7886982402	14239.1752748275
Organismal Systems	523.260523909263	510.797576381562	486.108494444667	532.243885256029	504.620819305394	488.969526609376	507.790434277258	469.104740826216	513.259867095318	477.373032928538	498.049487029437	527.389757311023	503.856194933095	501.730905570601	511.385253291288	496.360132167328	441.012404278403	499.936204060133
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

