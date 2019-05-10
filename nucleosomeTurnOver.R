library("ggplot2")
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst")
library("ggplot2")
#Figure1
setwd("./Figures/Fig1")
data<-read.delim(file="H3.H2BGFP.1kb.sum.tsv",header=F)
data<-na.omit(data[,c(1,2)])
cor(data[,1],data[,2],method="spearman")# 0.818701
cor(data[,1],data[,2])# 0.7162691
ggplot(data,aes(x=V1, y=V2))+stat_bin_2d(bins=500)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,2000),na.value="red")+labs(x="H3 signal",y="week0 H2BGFP signal")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,5))+scale_x_continuous(limits = c(0,5))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))
ggplot(data,aes(x=V1, y=V2))+stat_density_2d(geom="tile",aes(fill=..density..),contour=F,n = 100)+scale_fill_gradientn(colours = c("white","blue","yellow","red"))+
  scale_y_continuous(limits=c(0,5))+scale_x_continuous(limits = c(0,5))+geom_abline(slope=45,intercept = 0,linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
#Figure 2
setwd("./Figures/Fig2")
data<-read.delim(file="plus1e-4/NTR.summary.tsv",header=F)
data<-data[,7]
data[data>1]=1
hist<-hist(data[data>0],breaks=100,plot=F)
hist$counts<-hist$counts/sum(hist$counts)
plot(hist, freq=TRUE, ylab="Frequency",xlab="NTR",col="gray",main="")
data<-read.delim(file="plus1e-4/NTR.p0.05.txt",header=F)
hist<-hist(data$V1,breaks=100,plot=F)
hist$counts<-hist$counts/sum(hist$counts)
plot(hist, freq=TRUE, ylab="Frequency",xlab="NTR",col="gray",main="")
##Pie chart
library("RColorBrewer")
data<-read.delim(file="plus1e-4/top10k.p0.05.Anno.summary.tsv",header=F)
percent<-paste(round(data[,2]/sum(data[,2])*100,digits = 2),"%",sep="")
pie(data[,2],labels=paste(data[,1]," (",percent,")",sep=""),col=brewer.pal(n = 12, name = "Paired")[c(2,1,3,6,7,10,12)])
random<-read.delim(file="plus1e-4/random1k.Anno.summary.tsv",header=F)
fisher.test(matrix(c(2736,7264,710,9290),byrow = F,nrow = 2))#p-value < 2.2e-16
fisher.test(matrix(c(1004,8996,470,9530),byrow = F,nrow = 2))#p-value < 2.2e-16
data<-cbind(data$V2,random$V2)
colnames(data)<-c("Top10k","random10k")
rownames(data)<-random$V1
data[,1]<-data[,1]/sum(data[,1])
data[,2]<-data[,2]/sum(data[,2])
barplot(data,horiz=T,col=brewer.pal(n = 12, name = "Paired")[c(2,1,3,6,7,10,12)],legend=rownames(data),xlab="Percentage")

##Figure 2D
data<-read.delim(file="TSS.fl500.H2BGFP.fpkmSort.tsv",header=F)
library("plotrix")
logData<-log2(data[,c(2:7)]+0.0001)
colnames(logData)<-c("0w","1w","2w","4w","6w","8w")
logData<-as.matrix(logData)
q<-quantile(as.numeric(logData),probs=seq(0,1,0.1))
logData[logData<q[2]]=q[2]
logData<-rescale(logData,newrange = c(-1,1))
pheatmap(logData,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks=seq(-1,0,length.out = 6),color = colorRampPalette(c("darkblue","white","red"))(5))
data<-melt(logData)
ggplot(data, aes(Var2, Var1,fill=value)) + geom_tile()+scale_fill_gradientn(colors=c("blue3","white","red"))+
  labs(x="",y="")+theme_grey(base_size = -9)
fpkm<-as.matrix(log2(data$V8+1))
pheatmap(fpkm,cluster_rows = F,cluster_cols = F,show_rownames = F,breaks=seq(0,8,length.out = 51),color = colorRampPalette(rev(rainbow(7)))(50))
##Figure 2E
data<-read.delim(file="TSS.fl500.H2BGFP.fpkmSort.NTR.tsv",header=F)
bin<-c(rep(4,6242),rep(3,6242),rep(2,6242),rep(1,6243))
subset<-as.data.frame(cbind(bin,data$V9))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)
ggplot(subset,aes(x=FPKM, y=NTR, fill=FPKM))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2, lwd=1,fill=NA)+labs(y="NTR",x="Gene Expression(FPKM)")+
  ylim(-0.2,0.8)+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

wilcox.test(subset[subset$FPKM==1,2],subset[subset$FPKM==2,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==2,2],subset[subset$FPKM==3,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==3,2],subset[subset$FPKM==4,2],alternative = "l")#p-value < 2.2e-16
##Figure 2F
data<-read.delim(file="HM.lines.NTR.tsv",header=T)
step<-seq(-5000,5000,20)
colors<-rainbow(7)
plot(step,data[,1],ylim=c(-0.05,0.6),xlim=c(-5000,5000),xaxt="n",lwd=2,type="l",col=colors[1],xlab="Peak center relative position(bp)",ylab="NTR")
axis(1,at=seq(-5000,5000,1000),labels=seq(-5000,5000,1000))
for(i in 2:ncol(data)){
  lines(step,data[,i],col=colors[i],lwd=2)
}
legend("topright",legend=colnames(data),col=colors,lty=1,lwd=2)
data<-as.matrix(data)
data<-data-t(replicate(501,data[1,]))
plot(step,data[,1],ylim=c(-0.1,0.4),xlim=c(-5000,5000),xaxt="n",lwd=2,type="l",col=colors[1],xlab="Peak center relative position(bp)",ylab="NTR")
H3K4me1<-read.delim(file="H3K4me1.NTR.tsv",header=F)
H3K4me3<-read.delim(file="H3K4me3.NTR.tsv",header=F)
H3K27ac<-read.delim(file="enhancers.H2BGFP.NTR.tsv",header=F)
H3K27me3<-read.delim(file="H3K27me3.peaks.NTR.tsv",header=F)
H3K9ac<-read.delim(file="H3K9ac.NTR.tsv",header=F)
H3K36me3<-read.delim(file="H3K36me3.NTR.tsv",header=F)
H3K79me2<-read.delim(file="H3K79me2.NTR.tsv",header=F)
subset<-list(H3K4me3$V10,H3K4me1$V10,H3K9ac$V10,H3K27ac$V10,H3K27me3$V10,H3K36me3$V10,H3K79me2$V10)
names(subset)<-c("H3K4me1","H3K4me3","H3K9ac","H3K27ac","H3K27me3","H3K36me3","H3K79me2")
subset<-melt(subset)
subset$L1<-factor(subset$L1,levels=c("H3K4me1","H3K4me3","H3K9ac","H3K27ac","H3K27me3","H3K36me3","H3K79me2"))
ggplot(subset,aes(x=L1, y=value, fill=L1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.4, lwd=1,fill=NA,col=rainbow(7))+labs(y="NTR",x="")+
  ylim(-0.3,0.6)+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
interg<-read.delim(file="mm10.refGene.intergenic.1k.NTR.tsv",header=F)
subset<-list(H3K4me3$V10,H3K4me1$V10,H3K27ac$V10,H3K27me3$V10,interg$V10)
names(subset)<-c("H3K4me1","H3K4me3","H3K27ac","H3K27me3","intergenic")
subset<-melt(subset)
subset$L1<-factor(subset$L1,levels=c("H3K4me1","H3K4me3","H3K27ac","H3K27me3","intergenic"))
ggplot(subset,aes(x=L1, y=value, fill=L1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2, lwd=1,fill=NA,outlier.shape = NA,col=c("#FF0000FF","#FFDB00FF","#00FF92FF","#0092FFFF","black"))+labs(y="NTR",x="")+
  ylim(-0.3,0.6)+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
ggplot(subset,aes(x=L1, y=value, fill=L1))+geom_boxplot(fill=NA,outlier.shape = NA,col=c("#FF0000FF","#FFDB00FF","#00FF92FF","#0092FFFF","black"))+labs(y="NTR",x="")+
  ylim(-0.3,1.5)+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
#pvalues: H3K4me1/H3K4me3/H3K27ac/H3K27me3 vs. inter < 2.2e-16; H3K4me1/H3K4me3/H3K27ac vs. H3K27me3 < 2.2e-16; 
##Figure 2G
data<-read.delim(file="enhancers.TSSfl2k.closest.one-one.bed",header=F)
rst<-summary(unique(data[,c(15,16)])[,2])
bin1<-cbind(rep(1,nrow(data[data$V16<=rst[2],])),data[data$V16<=rst[2],10])
bin2<-cbind(rep(2,nrow(data[data$V16>rst[2] & data$V16<=rst[3],])),data[data$V16>rst[2] & data$V16<=rst[3],10])
bin3<-cbind(rep(3,nrow(data[data$V16>rst[3] & data$V16<=rst[5],])),data[data$V16>rst[3] & data$V16<=rst[5],10])
bin4<-cbind(rep(4,nrow(data[data$V16>rst[5],])),data[data$V16>rst[5],10])
subset<-as.data.frame(rbind(bin1,bin2,bin3,bin4))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)
ggplot(subset,aes(x=FPKM, y=NTR, fill=FPKM))+geom_violin(trim=F)+geom_boxplot(width=0.2, fill="white")+scale_fill_brewer(
  palette="Blues")+ylim(-0.3,0.8)+labs(y="NTR",x="Gene expression")+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$FPKM==1,2],subset[subset$FPKM==2,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==2,2],subset[subset$FPKM==3,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==3,2],subset[subset$FPKM==4,2],alternative = "l")#p-value < 2.2e-16
#eRNAs
data<-read.delim(file="enhancers.fpkm.NTR.tsv",header=F)
bin<-c(rep(1,6833),rep(2,6833),rep(3,6833),rep(4,6833))
subset<-as.data.frame(cbind(bin,data$V5))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)
summary(data$V4)
wilcox.test(data[data$V4==0,5],data[data$V4>1,5])
data<-read.delim(file="enhancers.fpkm.NTR.rmGeneRegion.tsv",header=F)
data<-data[order(data$V4),]
bin<-c(rep(1,1546),rep(2,1546),rep(3,1546),rep(4,1546))
boxplot(data[data$V4>0.2 & data$V4<=1,5],data[data$V4>1,5],ylim=c(-0.2,0.5),outline = F,names=c("Low","High"),xlab="eRNA expression",ylab="NTR",col="gray")
wilcox.test(data[data$V4>0.2 & data$V4<1,5],data[data$V4>=1,5])

data<-read.delim(file="eRNA.regions.RPKM.NTR.tsv",header=F)
data<-data[order(data$V4,decreasing = F),]
data<-data[data$V4>0.1,]
bin1<-cbind(rep(1,nrow(data[data$V4<=0.2,])),data[data$V4<=0.2,5])
bin2<-cbind(rep(2,nrow(data[data$V4>0.2 & data$V4<=1,])),data[data$V4>0.2 & data$V4<=1,5])
bin3<-cbind(rep(3,nrow(data[data$V4>1 & data$V4<=10,])),data[data$V4>1 & data$V4<=10,5])
bin4<-cbind(rep(4,nrow(data[data$V4>=10,])),data[data$V4>=10,5])
subset<-as.data.frame(rbind(bin1,bin2,bin3,bin4))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)#p value: 0.04,0.04,0.01


data<-read.delim()
data<-read.delim(file="enhancers.eRNA.NTR.forAggLines.tsv",header=F)
FPKM<-data[,c(1:5)]
data<-data[,c(6:ncol(data))]
plot(smooth.spline(colMeans(data[FPKM$V4>1,]),df = 100),type="l")
lines(colMeans(data[FPKM$V4<=0.1,]),col="red")

data<-read.delim(file="enahncersMACS14/enhancers.center1k.TSSfl2k.closest.one-one.bed",header=F)
bin<-c(rep(4,5687),rep(3,5687),rep(2,5687),rep(1,5689))
subset<-as.data.frame(cbind(bin,data$V4))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)

data<-read.delim(file="enhancersShareRegion/enhancers.TSSfl2k.closest.one-one.bed",header=F)
data<-read.delim(file="enhancersCenter/enhancers.center1k.TSSfl2k.closest.bed",header=F)
data<-read.delim(file="enhancersCenter/enhancers.center0.5k.TSSfl2k.closest.one-one.bed",header=F)
##Figure 2H & I
promoter<-read.delim(file="promoters.pol2.status.NTR.tsv",header=F)
enhancer<-read.delim(file="enhancers.pol2.status.NTR.final.tsv",header=F)
data<-matrix(nrow=4,ncol=3)
data[,1]<-c("promoters","promoters","enhancers","enhancers")
data[,2]<-c("wt","wo","wt","wo")
data[,3]<-c(nrow(promoter[promoter$V7=="wt",])/nrow(promoter),nrow(promoter[promoter$V7=="wo",])/nrow(promoter),
            nrow(enhancer[enhancer$V4=="wt",])/nrow(enhancer),nrow(enhancer[enhancer$V4=="wo",])/nrow(enhancer))
data<-as.data.frame(data,stringsAsFactors=F)
colnames(data)=c("class","pol2","Fraction")
data$class<-factor(data$class,levels=c("promoters","enhancers"))
data$pol2<-factor(data$pol2,levels=c("wt","wo"))
data$Fraction<-round(as.numeric(data$Fraction),digits = 2)
ggplot(data,aes(x=class,y=Fraction,fill=pol2))+geom_bar(stat="identity")+labs(y="Fraction")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
#promoters   wt 0.39
#promoters   wo 0.61
#enhancers   wt  0.26
#enhancers   wo  0.74
sub1<-cbind(rep("promoters",nrow(promoter)),promoter[,c(7,8)])
sub2<-cbind(rep("enhancers",nrow(enhancer)),enhancer[,c(4,5)])
colnames(sub1)<-c("class","pol2","NTR")
colnames(sub2)<-c("class","pol2","NTR")
subset<-rbind(sub1,sub2)
subset$pol2<-factor(subset$pol2,levels=c("wt","wo"))
ggplot(subset,aes(x=class, y=NTR, fill=pol2))+geom_violin(trim=F,position=position_dodge(0.9))+geom_boxplot(width=0.4,position=position_dodge(0.9))+ylim(0,0.7)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$class=="promoters" & subset$pol2=="wt",3],subset[subset$class=="promoters" & subset$pol2=="wo",3],alternative = "g")#p-value < 2.2e-16
wilcox.test(subset[subset$class=="enhancers" & subset$pol2=="wt",3],subset[subset$class=="enhancers" & subset$pol2=="wo",3],alternative = "g")#p-value < 2.2e-16
###Only enhancers
data<-as.data.frame(enhancer[,c(4,5)])
ggplot(data,aes(x=V4, y=V5, fill=V4))+geom_boxplot(width=0.4,position=position_dodge(0.9))+ylim(-0.5,1.0)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

## Figure 3
setwd("../fig3/")
promoter<-read.delim(file="promoters.TFs.status.NTR.tsv",header=F)
enhancer<-read.delim(file="enhancers.TFs.status.NTR.tsv",header=F)
data<-matrix(nrow=4,ncol=3)
data[,1]<-c("promoters","promoters","enhancers","enhancers")
data[,2]<-c("wt","wo","wt","wo")
data[,3]<-c(nrow(promoter[promoter$V7=="wt",])/nrow(promoter),nrow(promoter[promoter$V7=="wo",])/nrow(promoter),
            nrow(enhancer[enhancer$V4=="wt",])/nrow(enhancer),nrow(enhancer[enhancer$V4=="wo",])/nrow(enhancer))
data<-as.data.frame(data,stringsAsFactors=F)
colnames(data)=c("class","TFs","Fraction")
data$class<-factor(data$class,levels=c("promoters","enhancers"))
data$TFs<-factor(data$TFs,levels=c("wt","wo"))
data$Fraction<-round(as.numeric(data$Fraction),digits = 2)
ggplot(data,aes(x=class,y=Fraction,fill=TFs))+geom_bar(stat="identity")+labs(y="Fraction")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
#promoters   wt 0.18
# promoters   wo 0.82
# enhancers   wt 0.38
# enhancers   wo 0.62
sub1<-cbind(rep("promoters",nrow(promoter)),promoter[,c(7,8)])
sub2<-cbind(rep("enhancers",nrow(enhancer)),enhancer[,c(4,5)])
enhancer<-enhancer[enhancer$V3-enhancer$V2>20,]
colnames(sub1)<-c("class","TFs","NTR")
colnames(sub2)<-c("class","TFs","NTR")
subset<-rbind(sub1,sub2)
subset$TFs<-factor(subset$TFs,levels=c("wt","wo"))
ggplot(subset,aes(x=class, y=NTR, fill=TFs))+geom_violin(trim=F,position=position_dodge(0.9))+geom_boxplot(width=0.2,position=position_dodge(0.9))+ylim(0,0.7)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
ggplot(subset,aes(x=class, y=NTR, fill=TFs))+geom_boxplot(width=0.6,position=position_dodge(0.9))+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  +   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+ylim(-0.5,1)
wilcox.test(subset[subset$class=="promoters" & subset$TFs=="wt",3],subset[subset$class=="promoters" & subset$TFs=="wo",3],alternative = "g")#p-value < 2.2e-16
wilcox.test(subset[subset$class=="enhancers" & subset$TFs=="wt",3],subset[subset$class=="enhancers" & subset$TFs=="wo",3],alternative = "g")#p-value < 2.2e-16
###Only enhancers
data<-as.data.frame(enhancer[,c(4,5)])
ggplot(data,aes(x=V4, y=V5, fill=V4))+geom_boxplot(width=0.4,position=position_dodge(0.9))+ylim(-0.5,1.0)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))


promoter<-read.delim(file="promoters.TFnum.tsv",header=F)
boxplot(V8~V9, data = promoter[promoter$V9>0,],outline=F,xlab="TF number",ylab="NTR")
text(c(1,2,3,4),c(0.2,0.2,0.2,0.2),labels = c(2434,325,37,3))
promoter<-read.delim(file="allPromoters.TFnum.tsv",header=F)
boxplot(V8~V9, data = promoter[promoter$V9>0,],outline=F,xlab="TF number",ylab="NTR")
text(c(1.5,2.5,3.5),c(0.5,0.5,0.5),labels = c("p<2.2e-16","p=0.009","p=0.78"))

enhancer<-read.delim(file="enhancers.TFnum.NTR.finla.tsv",header=F)
enhancer<-enhancer[enhancer$V3-enhancer$V2>20,]
enhancer$V4<-enhancer$V4-1
subset<-as.data.frame(enhancer[,c(4,5)])
colnames(subset)<-c("status","NTR")
subset$status<-factor(subset$status)
subset[subset$status==4,1]=3
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_boxplot(lwd=1,fill=rgb(0,c(4,3,1.2)/4,0))+labs(x="TFs binding number",y="NTR")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  stat_summary(fun.y=mean, geom="point", size=2,color="white")
wilcox.test(subset[subset$status=="1",2],subset[subset$status=="2",2],alternative = "l")#1.361e-10
wilcox.test(subset[subset$status=="2",2],subset[subset$status=="3",2],alternative = "l")#0.03
#signal
data<-read.delim(file="enhancers.TFnum.allTFsignal.NTR.tsv",header=F)
data<-data[data$V3-data$V2>20,]
data[data$V4<0,4]=0
data[data$V5<0,5]=0
data[data$V6<0,6]=0
data[data$V7<0,7]=0
data<-data[order(rowMeans(data[,c(4:7)]),decreasing = F),]
bin<-c(rep(1,6492),rep(2,4359),rep(3,4109))#[0,10],(10,20),[20,)]
subset<-as.data.frame(cbind(bin,data$V9))
colnames(subset)<-c("Signals","NTR")
subset$Signals<-factor(subset$Signals)
ggplot(subset,aes(x=Signals, y=NTR, fill=Signals))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-1,1.6)+labs(y="NTR",x="TF signals")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$Signals=="1",2],subset[subset$Signals=="2",2],alternative = "l")#0.02;4.419e-07
#bin<-c(rep(1,nrow(data[data$V5<=5,])),rep(2,nrow(data[data$V5<=10 & data$V5>5,])),rep(3,nrow(data[data$V5<=20 & data$V5>10,])),rep(4,nrow(data[data$V5>20,])))
data<-as.data.frame(cbind(data$V8,rowMeans(data[,c(4:7)]),data$V9))
data[data$V2>60,2]<-60
scatter3D(data$V1,data$V2,data$V3,pch=19,xlim=c(2,6),ylim=c(-1,60),bty="u",col.panel =NA,theta = 15,phi=30,xlab="TFs binding number",ylab="Total TF signals",zlab="NTR",ticktype = "detailed")
data2<-data
data2$V2=floor(data2$V2/2.5)*2.5
scatter3D(data2$V1,data2$V2,data2$V3,pch=19,xlim=c(2,6),bty="u",col.panel =NA,theta = 15,phi=30,xlab="TFs binding number",ylab="Total TF signals",zlab="NTR",ticktype = "detailed")

data<-read.delim(file="Gata4.Nkx2-5.Tbx20.Tbx3.agg.txt",header=F)
step=seq(-3000,3000,by=20)
colnames(data)<-c("Gata4","Nkx2-5","Tbx20","Tbx3")
data<-data[c(100:400),]
colors=rainbow(n=4)
plot(step,data[,1],col=colors[1],type="l",xlab="Peak center relative position(bp)",ylab="HTR",xlim=c(-3000,3000),ylim=c(0,0.6),lwd=2)
lines(step,data[,2],col=colors[2],lwd=2)
lines(step,data[,3],col=colors[3],lwd=2)
lines(step,data[,4],col=colors[4],lwd=2)
legend("topright",legend=colnames(data),col=colors,lty=1,lwd=2)
#All enhancers 20180924
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/Figures/Fig3/allEnhancers")
enhancer<-read.delim(file="enhancers.TFs.status.NTR.tsv",header=F)
data<-matrix(nrow=2,ncol=3)
data[,1]<-c("enhancers","enhancers")
data[,2]<-c("wt","wo")
data[,3]<-c(0.5222,0.4778)#wt:14274,wo:13058
data<-as.data.frame(data,stringsAsFactors=F)
colnames(data)=c("class","TFs","Fraction")
data$TFs<-factor(data$TFs,levels=c("wt","wo"))
data$Fraction<-round(as.numeric(data$Fraction),digits = 2)
ggplot(data,aes(x=class,y=Fraction,fill=TFs))+geom_bar(stat="identity")+labs(y="Fraction")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
data<-as.data.frame(enhancer[,c(4,5)])
ggplot(data,aes(x=V4, y=V5, fill=V4))+geom_boxplot(width=0.4,position=position_dodge(0.9),outlier.shape = NA)+ylim(-0.5,1.0)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

enhancer<-read.delim(file="enhancers.TFnum.NTR.final.tsv",header=F)
enhancer<-enhancer[enhancer$V3-enhancer$V2>20,]
enhancer$V4<-enhancer$V4-1
subset<-as.data.frame(enhancer[,c(4,5)])
colnames(subset)<-c("status","NTR")
subset$status<-factor(subset$status)
subset[subset$status==4,1]=3
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_boxplot(lwd=1,fill="gray")+labs(x="TFs binding number",y="NTR")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  stat_summary(fun.y=mean, geom="point", size=2,color="white")
wilcox.test(subset[subset$status=="1",2],subset[subset$status=="2",2],alternative = "l")#<2.2e-16
wilcox.test(subset[subset$status=="2",2],subset[subset$status=="3",2],alternative = "l")#0.32

data<-read.delim(file="enhancers.TFnum.allTFsignal.NTR.tsv",header=F)
data<-data[data$V3-data$V2>20,]
data<-na.omit(data)
data[data$V4<0,4]=0
data[data$V5<0,5]=0
data[data$V6<0,6]=0
data[data$V7<0,7]=0
data<-data[order(rowMeans(data[,c(4:7)]),decreasing = F),]
bin<-c(rep(1,12247),rep(2,13658),rep(3,4300))#[0,10],(10,20),[20,)]
subset<-as.data.frame(cbind(bin,data$V9))
colnames(subset)<-c("Signals","NTR")
subset$Signals<-factor(subset$Signals)
ggplot(subset,aes(x=Signals, y=NTR, fill=Signals))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-1,1.6)+labs(y="NTR",x="TF signals")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$Signals=="1",2],subset[subset$Signals=="2",2],alternative = "l")#<2.2e-16; 2.822e-09

#Figure4_final
setwd("../Fig4/")
data<-read.delim(file="enhancers.NTR.H3K27ac.signal.tsv",header=F)
bin<-c(rep(1,6833),rep(2,6833),rep(3,6833),rep(4,6835))
subset<-as.data.frame(cbind(bin,data$V4))
colnames(subset)<-c("FPKM","NTR")
subset$FPKM<-factor(subset$FPKM)
ggplot(subset,aes(x=FPKM, y=NTR, fill=FPKM))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(0,0.6)+labs(y="NTR",x="Gene Expression(FPKM)")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$FPKM==1,2],subset[subset$FPKM==2,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==2,2],subset[subset$FPKM==3,2],alternative = "l")#p-value < 2.2e-16
wilcox.test(subset[subset$FPKM==3,2],subset[subset$FPKM==4,2],alternative = "l")#p-value < 2.2e-16
data<-read.delim(file="enhancers.H3K27ac.signal.H2BGFP0.subtract.NTR.tsv",header=F)
#Other histone modifications
data<-read.delim(file="../Fig2/signalsTest/H3K9ac.signals.NTR.tsv",header=F)
bin<-c(rep(1,11370),rep(2,11370),rep(3,11370),rep(4,11373))
subset<-as.data.frame(cbind(bin,data$V5))
data<-read.delim(file="../Fig2/signalsTest/H3K36me3.signals.NTR.tsv",header=F)
bin<-c(rep(1,15162),rep(2,15162),rep(3,15162),rep(4,15162))
data<-read.delim(file="../Fig2/signalsTest/H3K79me2.signals.NTR.tsv",header=F)
bin<-c(rep(1,8372),rep(2,8372),rep(3,8372),rep(4,8373))

setwd("./enhancerNarrowPeak/")
data<-read.delim(file="enhancers.NTR.H3K27ac.signal.tsv",header=F)
bin<-c(rep(1,9393),rep(2,9393),rep(3,9393),rep(4,9395))
data<-read.delim(file="enhancers.NTR.coBindNum.tsv",header=F)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(0,0.6)+labs(y="NTR",x="Regulators binding number")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
#EED and H3K27me3
h3k27me3<-read.delim(file="H3K27me3.unique.NTR.tsv",header=F)
share<-read.delim(file="H3K27me3.EED.shareRegion.NTR.tsv",header=F)
eed<-read.delim(file="EED.unique.NTR.tsv",header=F)
h3k27me3<-cbind(h3k27me3[,10],1)
share<-cbind(share[,10],2)
eed<-cbind(eed[,10],3)
subset<-as.data.frame(rbind(h3k27me3,share,eed))
share<-read.delim(file="H3K27me3.SUZ12.shareRegion.NTR.tsv",header=F)
h3k27me3<-read.delim(file="H3K27me3.SUZ12.intersect.bed",header=F)
suz12<-read.delim(file="SUZ12.unique.vsK27me3.NTR.tsv",header=F)
h3k27me3<-cbind(h3k27me3[h3k27me3$V6<0,4],1)
share<-cbind(share[,10],2)
suz12<-cbind(suz12[,10],3)
subset<-as.data.frame(rbind(h3k27me3,share,suz12))
##EED & H3K27me3 genes FPKM
h3k27me3<-read.delim(file="H3K27me3genesFPKM/H3K27me3.unique.vsEED.tsv",header=F)
share<-read.delim(file="H3K27me3genesFPKM/H3K27me3.EED.share.tsv",header=F)
eed<-read.delim(file="H3K27me3genesFPKM/EED.unique.tsv",header=F)
others<-read.delim(file="H3K27me3genesFPKM/noPeak.associated.genes.bed6",header=F)
all<-read.delim(file="../Fig2/TSS.fl500.H2BGFP.fpkmSort.tsv",header = F)
h3k27me3<-unique(h3k27me3[,c(4:9)])
share<-unique(share[,c(4:9)])
eed<-unique(eed[,c(4:9)])
h3k27me3<-cbind(h3k27me3[,5],1)
share<-cbind(share[,5],2)
eed<-cbind(eed[,5],3)
others<-cbind(others[,5],4)
all<-cbind(all[,8],5)
subset<-as.data.frame(rbind(h3k27me3,share,eed,others,all))
colnames(subset)<-c("FPKM","status")
subset$status<-factor(subset$status)
subset$FPKM<-log2(subset$FPKM+1)
ggplot(subset,aes(x=status, y=FPKM, fill=status))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.1,lwd=1,fill=NA,color=c("blue","red","purple","black","black"))+
     ylim(0,5)+labs(y="Genes FPKM",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_x_discrete(labels=c("H3K27me3+/EED-","H3K27me3+/EED+","H3K27me3-/EED+","H3K27me3-/EED-","all"))
boxplot(h3k27me3$V8,share$V8,eed$V8,outline = F,names = c("H3K27me3+/EED-","H3K27me3+/EED+","H3K27me3-/EED+"),ylab="Genes FPKM",ylim=c(0,20))
###peak regions transcription states
h3k27me3<-read.delim(file="H3K27me3genesFPKM/H3K27me3.unique.hete.repMean.RPKM.tsv",header=F)
eed<-read.delim(file="H3K27me3genesFPKM/EED.unique.hete.repMean.RPKM.tsv",header=F)
share<-read.delim(file="H3K27me3genesFPKM/shared.sorted.hete.repMean.RPKM.tsv",header=F)
h3k27me3<-cbind(h3k27me3[,4],1)
share<-cbind(share[,4],2)
eed<-cbind(eed[,4],3)
subset<-as.data.frame(rbind(h3k27me3,share,eed))

#final vesion
#EED & enhancer
enhancer<-read.delim(file="eed/EED.enhancers.intersect.V1.bed",header=F)
share<-read.delim(file="eed/EED.enhancers.shareRegion.V1.NTR.tsv",header=F)
regulator<-read.delim(file="eed/EED.unique.V1.NTR.tsv",header=F)
grid.arrange(H3K27me3,H3K27ac,ncol=2,nrow=1)

#HDAC1 & enhancer
enhancer<-read.delim(file="hdac1/HDAC1.enhancers.intersect.V1.bed",header=F)
share<-read.delim(file="hdac1/HDAC1.enhancers.shareRegion.V1.NTR.tsv",header=F)
regulator<-read.delim(file="hdac1/HDAC1.unique.V1.NTR.tsv",header=F)

enhancer<-read.delim(file="hdac2/HDAC2.enhancers.intersect.V1.bed",header=F)
share<-read.delim(file="hdac2/HDAC2.enhancers.shareRegion.V1.NTR.tsv",header=F)
regulator<-read.delim(file="hdac2/HDAC2.unique.V1.NTR.tsv",header=F)

enhancer<-read.delim(file="p300/p300.enhancers.intersect.V0.bed",header=F)
share<-read.delim(file="p300/p300.enhancers.shareRegion.V0.NTR.tsv",header=F)
regulator<-read.delim(file="p300/p300.unique.V0.NTR.tsv",header=F)

enhancer<-read.delim(file="suz12/suz12.enhancers.intersect.V2.bed",header=F)
share<-read.delim(file="suz12/suz12.enhancers.shareRegion.V2.NTR.tsv",header=F)
regulator<-read.delim(file="suz12/suz12.unique.V2.NTR.tsv",header=F)

enhancer<-cbind(enhancer[enhancer$V6<0,4],1)
share<-cbind(share[,10],2)
regulator<-cbind(regulator[,10],3)
subset<-as.data.frame(rbind(enhancer,share,regulator))
colnames(subset)<-c("NTR","status")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA,color=c("blue","red","purple"))+
  ylim(-0.2,0.6)+labs(y="NTR",x="SUZ12")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$status=="1",1],subset[subset$status=="2",1],alternative = "l")
wilcox.test(subset[subset$status=="3",1],subset[subset$status=="2",1],alternative = "l")
###10181016
enhancer<-cbind(enhancer[enhancer$V6<0,4],1)
share<-cbind(share[,10],2)
subset<-as.data.frame(rbind(enhancer,share))
colnames(subset)<-c("NTR","status")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_boxplot(lwd=1,fill=NA,color=c("blue","red"))+
  ylim(-0.2,1.5)+labs(y="NTR",x="")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
h3k27me3<-read.delim(file="H3K27me3.unique.NTR.tsv",header=F)
share<-read.delim(file="H3K27me3.EED.shareRegion.NTR.tsv",header=F)
h3k27me3<-cbind(h3k27me3[,10],1)
share<-cbind(share[,10],2)
subset<-as.data.frame(rbind(h3k27me3,share))

data<-read.delim(file="enhancers.coBind.NTR.tsv",header=F)
data$V4<-data$V4-1
subset<-as.data.frame(cbind(data[,5],data[,4]))
colnames(subset)<-c("NTR","status")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_boxplot(lwd=1,fill=rgb(0,c(4,3,2,1.2)/4,0))+labs(x="Regulators bindind number",y="NTR")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
stat_summary(fun.y=mean, geom="point", size=2,color="white")
###20181016
data<-read.delim(file="enhancers.coBind.NTR.tsv",header=F)
data$V4<-data$V4-1
data[data$V4>=3,4]=3
data<-data[data$V4>0,]
subset<-as.data.frame(cbind(data[,5],data[,4]),stringsAsFactors = F)
colnames(subset)<-c("NTR","status")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-1,2)+labs(y="NTR",x="# of co-occupied modifiers")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
data<-read.delim(file="H3K27me3.coBind.NTR.tsv",header = F)
data$V4<-data$V4-1
subset<-as.data.frame(cbind(data[,5],data[,4]),stringsAsFactors = F)
colnames(subset)<-c("NTR","status")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-1,2)+labs(y="NTR",x="# of co-occupied modifiers")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

boxplot(NTR~status,data=subset,outline=F,col=rgb(0,c(4,3,2,1)/4,0),xlab="Regulators binding number",ylab="NTR")
##random background
me3<-read.delim(file="H3K27me3.unique.NTR.tsv",header=F)
me3_random<-read.delim(file="H3K27me3.EEDnegative.random.NTR.tsv",header=F)
ac<-read.delim(file="eed/EED.enhancers.intersect.V1.bed",header=F)
ac<-ac[ac$V6<0,]
ac_random<-read.delim(file="eed/H3K27ac.EEDnegative.random.NTR.tsv",header = F)
vioplot(me3$V10,me3_random$V10,ac$V4,ac_random$V10)
data<-rbind(cbind("H3K27me3",me3$V10),cbind("random",na.omit(me3_random$V10)))
data<-as.data.frame(data,stringsAsFactors =F)
data$V2<-as.numeric(data$V2)
data$V1<-as.factor(data$V1)
ggplot(data,aes(x=V1, y=V2, fill=V1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-0.2,0.6)+labs(y="HTR",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
data<-rbind(cbind("H3K27ac",ac$V4),cbind("random",na.omit(ac_random$V10)))

#Figure5 
setwd("../Fig5/")
common<-read.delim(file="WT.EEDko.common.EEDko-WT.signals.tsv",header=F)
boxplot(common$V4,common$V5,outline=F,names = c("KO","WT"),ylab="H3K27ac signals")
data<-read.delim(file="1kbIntervals.H3K27ac.EEDko.WT.signals.tsv",header=F)
data<-data[data$V4>0 & data$V5>0,]
boxplot(data$V4,data$V5,outline=F,names = c("KO","WT"),ylab="H3K27ac signals")
common<-read.delim(file="commonPeaks.signals.danpos.tsv",header=F)
#Fig5B
setwd("./H3K27me3")
library("gridExtra")
unchange<-read.delim(file="EEDko.unchange.wt.ko.NTR.bed3+",header=F)
loss<-read.delim(file="EEDko.loss.wt.ko.NTR.bed3+",header=F)
boxplot(unchange$V4,unchange$V5,loss$V4,loss$V5,outline = F,ylab="NTR",names=c("unchange_wt","unchange_ko","loss_wt","loss_ko"),at=c(1,1.9,3,3.9))
text(c(1.5,3.5),c(0.4,0.4),c("***","***"))
unchange<-read.delim(file="EEDko.unchange.wt.ko.NTR.agg.txt",header = F)
loss<-read.delim(file="EEDko.loss.wt.ko.NTR.agg.txt",header=F)
plot(unchange$V1,unchange$V2,type="l",col="blue",xlab = "Peak center relative position",ylab="NTR",ylim=c(0,0.2))
lines(unchange$V1,unchange$V3,col="red")
legend("topright",c("EEDhete","EEDko"),col=c("blue","red"),lty=1)
plot(loss$V1,loss$V2,type="l",col="blue",xlab = "Peak center relative position",ylab="NTR",ylim=c(0,0.2))
lines(loss$V1,loss$V3,col="red")
legend("topright",c("EEDhete","EEDko"),col=c("blue","red"),lty=1)
loss<-read.delim(file="EEDko.loss.ko.wt.NTR.signalDiff.tsv",header=F)
unchange<-read.delim(file="EEDko.unchange.ko.wt.NTR.signalDiff.tsv",header=F)
subset<-rbind(loss[,c(4,5)],unchange[,c(4,5)])
ggplot(subset,aes(x=V5, y=V4))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="H3K27me3 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-20,5))
subset<-read.delim(file="../H3K27ac/peaks.ko.wt.signal.NTR.diff.tsv",header=F)
ggplot(subset,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red",limits=c(0,200))+labs(x="H3K27ac signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-20,10))
subset<-read.delim(file="../H3K27ac/peaks.ko.wt.signalFC.NTR.diff.tsv",header=F)
ggplot(subset,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,200))+labs(x="H3K27ac signals(CKO/WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-1,10))
subset<-read.delim(file="../H3K27me3/peaks.ko.wt.signalFC.NTR.diff.tsv",header=F)
scatterPlot<-ggplot(subset,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,200))+labs(x="H3K27me3 signals(CKO/WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-1,10))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))
densityx<-ggplot(subset, aes(x=V4)) +geom_density(alpha=.5,fill="gray")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_x_continuous(limits = c(-1,10))+labs(x="H3K27me3 signals(CKO/WT)")
densityy<-ggplot(subset, aes(x=V5)) +geom_density(alpha=.5,fill="gray")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_x_continuous(limits = c(-2,2))+labs(x="NTR(CKO-WT)")+coord_flip()
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank()
  )
grid.arrange(densityx, blankPlot,scatterPlot, densityy, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))


x<-matrix(c(20079,28646,8697,13701),nrow=2)
rownames(x)<-c("increase","decrease")
colnames(x)<-c("H3K27ac","H3K27me3")
barplot(x,ylab="Peak number",col=c(gray(0.9),gray(0.4)))
legend("topright",c("increase","decrease"),fill=c(gray(0.9),gray(0.4)))
fisher.test(x) #p-value = 1.833e-09
plot(ecdf(subset1$V5),xlim=c(-1,1),verticals = T, do.points=F,xlab="NTR(CKO-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(subset2$V5),col="red",verticals = T, do.points=F)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("black","red"), lty=1)
ks.test(subset1$V5,subset2$V5,alternative = "g")#p-value < 2.2e-16
text(-0.7,0.6,"KS test:p<2.2e-16")
###
subset1<-read.delim(file="./H3K27ac/peaks.ko.wt.repSignal.NTR.diff.tsv",header=F)
subset1[subset1$V4<0,4]=0
subset1[subset1$V5<0,5]=0
subset1[subset1$V6<0,6]=0
subset1[subset1$V7<0,7]=0
subset1<-as.data.frame(cbind(rowMeans(subset1[,c(4,5)])/rowMeans(subset1[,c(6,7)]),subset1$V8))
subset1<-subset1[subset1$V1>1.5,]
#subset1<-subset1[subset1$V1>1,]
H3K27ac<-ggplot(subset1,aes(x=V1, y=V2))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),na.value="red",limits=c(0,20))+labs(x="H3K27ac signals(CKO/WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(1.5,10),breaks = seq(0.5,10,1))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))+geom_abline(slope=0,intercept = 0,col="red",lwd=1)
subset2<-read.delim(file="./H3K27me3/peaks.ko.wt.repSignal.NTR.diff.tsv",header=F)
subset2[subset2$V4<0,4]=0
subset2[subset2$V5<0,5]=0
subset2[subset2$V6<0,6]=0
subset2[subset2$V7<0,7]=0
subset2[subset2$V8<0,8]=0
subset2<-as.data.frame(cbind(rowMeans(subset2[,c(4,5)])/rowMeans(subset2[,c(6,7,8)]),subset2$V9))
subset2<-subset2[subset2$V1<0.5,]
#subset2<-subset2[subset2$V1<1,]
H3K27me3<-ggplot(subset2,aes(x=V1, y=V2))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),na.value ="red",limits=c(0,20))+labs(x="H3K27me3 signals(CKO/WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))+geom_abline(slope=0,intercept = 0,col="red",lwd=1)
grid.arrange(H3K27ac,H3K27me3,ncol=1,nrow=2)
#density plot
subset<-rbind(cbind(subset1$V2,"H3K27ac"),cbind(subset2$V2,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,color=V2))+geom_density()+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KO-WT)")
#barplot
x<-matrix(c(5076,2380,22173,18811),nrow=2)
rownames(x)<-c("decrease","increase")
colnames(x)<-c("H3K27ac","H3K27me3")
x[,1]=x[,1]/sum(x[,1])
x[,2]=x[,2]/sum(x[,2])
barplot(x,ylab="Fraction",col=c(gray(0.4),gray(0.9)))
legend("topright",c("increase","decrease"),fill=c(gray(0.9),gray(0.4)))
x<-round(c(x[1,1],x[2,1],x[1,2],x[2,2]),digits = 4)*100
percent=""
percent[1]=paste(x[1],"%",sep="")
percent[2]=paste(x[2],"%",sep="")
percent[3]=paste(x[3],"%",sep="")
percent[4]=paste(x[4],"%",sep="")
text(c(0.8,0.8,1.8,1.8),c(0.3,0.8,0.3,0.8),percent)
fisher.test(x) #p-value <2.2e-16
#ECDF
plot(ecdf(subset1$V2),xlim=c(-1,1),verticals = T, do.points=F,xlab="NTR(CKO-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(subset2$V2),col="red",verticals = T, do.points=F)
abline(h=0.5,col="gray",lty=2)
abline(v=0,col="gray",lty=2)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("black","red"), lty=1)
ks.test(subset1$V2,subset2$V2)#p-value < 2.2e-16
text(-0.7,0.6,"KS test:p<2.2e-16")
#NTR binned by H3K27me3 signals
subset2<-read.delim(file="./H3K27me3/WT.peaks.signals.NTR.tsv",header=F)
subset2<-subset2[order(subset2$V4,decreasing = F),]
bin<-c(rep(1,13147),rep(2,13147),rep(3,13147),rep(4,13150))
subset<-as.data.frame(cbind(bin,subset2$V5))
subset$bin<-as.factor(subset$bin)
ggplot(subset,aes(x=bin, y=V2, fill=bin))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-0.1,0.8)+labs(y="NTR",x="H3K27me3 signals")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$bin=="1",2],subset[subset$bin=="2",2],alternative = "l")
wilcox.test(subset[subset$bin=="2",2],subset[subset$bin=="3",2],alternative = "l")
wilcox.test(subset[subset$bin=="3",2],subset[subset$bin=="4",2],alternative = "l")
data<-read.delim(file="./H3K27me3/WT.peaks.signal.H2BGFP0.subtract.NTR.tsv",header=F)
bin<-c(rep(1,13147),rep(2,13147),rep(3,13147),rep(4,13150))
subset<-as.data.frame(cbind(bin,data$V7))
data<-read.delim(file="./H3k27me3/WT.merged.peaks.signals.NTR.tsv",header=F)
data<-data[order(data$V4,decreasing = F),]
bin<-c(rep(1,46546),rep(2,46546),rep(3,46546),rep(4,46549))
subset<-as.data.frame(cbind(bin,data$V5))
data<-read.delim(file="./H3K27me3/WT.merged.peaks.signals.NTR.FCincrease.tsv",header=F)

###
#nucleosome occupancy
data<-read.delim(file="./H3K27me3/peaks.wt.ko.NC.tsv",header=F)
agg<-read.delim(file="./H3K27me3/peaks.NC.agg.txt",header=F)
plot(agg$V1,(agg$V2+agg$V3)/2,type="l",col="red",xlab="Peak relative position",ylab="Normalized nucleosome occupancy")
lines(agg$V1,(agg$V4+agg$V5)/2,col="blue")
legend("topleft",legend = c("WT","EEDko"),col=c("red","blue"),lty=1)
agg<-read.delim(file="./H3K27ac/peaks.NC.agg.txt",header=F)
plot(agg$V1,(agg$V2+agg$V3)/2,type="l",col="red",ylim=c(0.25,0.35),xlab="Peak relative position",ylab="Normalized nucmeosome occupancy")
lines(agg$V1,(agg$V4+agg$V5)/2,col="blue")
legend("topright",legend = c("WT","EEDko"),col=c("red","blue"),lty=1)
plot(agg$V1,agg$V3,type="l",col="red",ylim=c(0.25,0.4),xlab="Peak relative position",ylab="Normalized nucmeosome occupancy")
lines(agg$V1,agg$V2,col="brown")
lines(agg$V1,agg$V4,col="blue")
lines(agg$V1,agg$V5,col="darkblue")
legend("topright",legend = c("WT-1","WT-2","EEDko-1","EEDko-2"),col=c("brown","red","blue","darkblue"),lty=1)

agg<-read.delim(file="./H3K27ac/peaks.ko.wt.signalFC.gt1.5.NCagg.txt",header=F)
agg<-read.delim(file="./H3K27me3/peaks.ko.wt.signalsFC.lt0.5.NCagg.txt",header=F)
par(mfrow=c(2,1))
data<-as.matrix(agg[,c(2:5)])
#Norm by subtract
rst<-data-t(replicate(10000,data[1,]))
plot(agg[,1],(rst[,1]+rst[,2])/2,type="l",col="red",ylim=c(-0.1,0.1),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
lines(agg[,1],(rst[,3]+rst[,4])/2,col="blue")
legend("topright",legend = c("WT","EEDko"),col=c("red","blue"),lty=1)
#Norm by whole region mean(Final version)
rst<-data/t(replicate(10000,colMeans(data)))
plot(agg[,1],(rst[,1]+rst[,2])/2,type="l",col="red",ylim=c(0.8,1.2),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
lines(agg[,1],(rst[,3]+rst[,4])/2,col="blue")
legend("topright",legend = c("WT","EEDko"),col=c("red","blue"),lty=1)
wilcox.test(((rst[,1]+rst[,2])/2)[c(4000:6000)],((rst[,3]+rst[,4])/2)[c(4000:6000)],paired = T)#p-value  < 2.2e-16 for H3K27ac & p - value p-value < 2.2e-16 for H3K27me3
plot(agg[,1],rst[,1],type="l",col="red",ylim=c(0.8,1.2),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
lines(agg[,1],rst[,2],col="brown")
lines(agg[,1],rst[,3],col="blue")
lines(agg[,1],rst[,4],col="blue")
legend("topright",legend = c("WT-1","WT-2","EEDko-1","EEDko-2"),col=c("brown","red","darkblue","blue"),lty=1)
t.test(c(mean(rst[c(4000:6000),1]),mean(rst[c(4000:6000),2])),c(mean(rst[c(4000:6000),3]),mean(rst[c(4000:6000),4])))#0.7876 for H3K27me3; 0.7739 for H3K27ac
#Norm by center 2000 mean
rst<-data/t(replicate(10000,colMeans(data[c(4000,6000),])))
plot(agg[,1],(rst[,1]+rst[,2])/2,type="l",col="red",ylim=c(0.8,1.2),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
#Norm by center 4000 mean
rst<-data/t(replicate(10000,colMeans(data[c(3000,7000),])))
#Normlized GC
agg<-read.delim(file="./H3K27ac/peaks.ko.wtFC0.5.GCnorm.NC.agg.txt",header=F)

data<-read.delim(file="./H3K27ac/peaks.ko.wt.signalFC.gt0.5.NC.tsv",header=F)
data<-cbind(rowMeans(data[,c(4,5)]),rowMeans(data[,c(6,7)]))
data<-as.data.frame(data)
ggplot(data,aes(x=V1, y=V2))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,150))+labs(x="Normalized nucleosome occupancy in WT",y="Normalzied nucleosome occupancy in EED KO")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits = c(0,1))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))+geom_abline(slope=1,intercept = 0,col="red",lwd=1)
#Figure6 
setwd("../Fig6")
unchange<-read.delim(file="NTR/H3K27acNarrowPeaks/banding_sham.unchange.H3K27ac.summary.tsv",header=F)
banding_induce<-read.delim(file="NTR/H3K27acNarrowPeaks/banding_increase.H3K27ac.summary.tsv",header=F)
banding_loss<-read.delim(file="NTR/H3K27acNarrowPeaks/sham_increase.H3K27ac.summary.tsv",header=F)
unchange<-cbind(rowMeans(unchange[,c(2,3)]),rowMeans(unchange[,c(4,5)]))
colnames(unchange)<-c("banding","sham")
banding_induce<-cbind(rowMeans(banding_induce[,c(2,3)]),rowMeans(banding_induce[,c(4,5)]))
colnames(banding_induce)<-c("banding","sham")
banding_loss<-cbind(rowMeans(banding_loss[,c(2,3)]),rowMeans(banding_loss[,c(4,5)]))
colnames(banding_loss)<-c("banding","sham")
unchange<-unchange*(-1)
banding_loss<-banding_loss*(-1)
banding_induce<-banding_induce*(-1)
unchange<-melt(unchange)[,c(2,3)]
unchange<-cbind(replicate(nrow(unchange),"unchange"),unchange)
banding_induce<-melt(banding_induce)[,c(2,3)]
banding_induce<-cbind(replicate(nrow(banding_induce),"banding_induce"),banding_induce)
banding_loss<-melt(banding_loss)[,c(2,3)]
banding_loss<-cbind(replicate(nrow(banding_loss),"banding_loss"),banding_loss)
colnames(unchange)<-c("class","condition","value")
colnames(banding_loss)<-c("class","condition","value")
colnames(banding_induce)<-c("class","condition","value")
data<-rbind(unchange,banding_loss,banding_induce)
ggplot(data, aes(x = class, y = value, fill = condition)) +geom_boxplot(outlier.shape = NA) + ylim(-0.5,0) +labs(x="",y="Relative nucleosome turnover")+theme(axis.text.x = element_text(angle = 90))

setwd("./H3K27acSignals/")
sham<-read.delim(file="sham.agg.txt",header=F)
banding<-read.delim(file="banding.agg.txt",header=F)
plot(banding$V1,banding$V2,type="l",col="red",lwd=2,ylab="H3K17ac ChIP-seq signals",xlab="Peak center relative position")
lines(sham$V1,sham$V2,col="blue",lwd=2)
legend("topright",legend = c("banding","sham"),col=c("red","blue"),lty=1,lwd=2)
##Other data processing
#MNase-seq 
setwd("E:/百度云同步盘/my PC D/projects/help/nucleosome/NewRst/MNase-seq")
wt_1<-read.delim(file="1_MNase-EEDhetoMF-lib1610/tagGCcontent.txt",header=T)
wt_2<-read.delim(file="1_MNase-EEDhetoMF-lib1612/tagGCcontent.txt",header=T)
ko_1<-read.delim(file="1_MNase-EEDkoMF-lib1611/tagGCcontent.txt",header=T)
ko_2<-read.delim(file="1_MNase-EEDkoMF-lib1613/tagGCcontent.txt",header=T)
plot(wt_1[,1],wt_1[,3],type="l",col="red",xlab="GC content of fragments",ylab="Normalized fraction",ylim=c(0,6),lwd=2)
lines(wt_2[,1],wt_2[,3],col="brown",lwd=2)
lines(ko_1[,1],ko_1[,3],col="blue",lwd=2)
lines(ko_2[,1],ko_2[,3],col="darkblue",lwd=2)
legend("topright",legend = c("WT-1","WT-2","EEDko-1","EEDko-2"),col=c("red","brown","blue","darkblue"),lty=1)

wt_1<-read.delim(file="1_MNase-EEDhetoMF-lib1610_correct/tagGCnormalization.txt",header=T)
wt_2<-read.delim(file="1_MNase-EEDhetoMF-lib1612_correct/tagGCnormalization.txt",header=T)
ko_1<-read.delim(file="1_MNase-EEDkoMF-lib1611_correct/tagGCnormalization.txt",header=T)
ko_2<-read.delim(file="1_MNase-EEDkoMF-lib1613_correct/tagGCnormalization.txt",header=T)

##Tesyt signals
setwd("E:/百度云同步盘/my PC D/projects/help/nucleosome/NewRst/Figures/Fig4/NTR_by_Signals")
files<-list.files(pattern="H3K27ac")
for(i in 1:4){
  data<-read.delim(file=files[i],header=F)
  #bin<-c(rep(1,6833),rep(2,6833),rep(3,6833),rep(4,6833))
  bin<-c(rep(1,13147),rep(2,13147),rep(3,13147),rep(4,13150))
  subset<-as.data.frame(cbind(bin,data$V5))
  colnames(subset)<-c("FPKM","NTR")
  subset$FPKM<-factor(subset$FPKM)
  ggplot(subset,aes(x=FPKM, y=NTR, fill=FPKM))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
    ylim(-0.2,0.6)+labs(y="NTR",x=files[i])+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
}
files=list.files(pattern = "H3K27me3")
#NTR & H3K27ac & H3K27me3
setwd("E:/百度云同步盘/my PC D/projects/help/nucleosome/NewRst/Figures/Fig4/basedNucleosomes")
data<-read.delim(file="nucleosome.H3K27ac.H3K27me3.NTR.tsv",header=F)
subset<-data[data$V7<0.05 & data$V8>0.8,c(4,5,6)]
subset<-subset[order(subset$V6,decreasing = F),]
bin<-c(rep(1,10997),rep(2,10997),rep(3,10997),rep(4,10999))
subset<-as.data.frame(cbind(bin,subset$V4))

data<-data[order(data$V7,decreasing = F),]
data1<-data[data$V4>=0,]
bin<-c(rep(1,43280),rep(2,43280),rep(3,43280),rep(4,43281))
subset<-as.data.frame(cbind(bin,data1$V6))
colnames(subset)<-c("signals","NTR")
subset$signals<-factor(subset$signals)
ggplot(subset,aes(x=signals, y=NTR, fill=signals))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-0.2,0.6)+labs(y="signals",x="NTR")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$signals==1,2],subset[subset$signals==2,2],alternative="l")
wilcox.test(subset[subset$signals==2,2],subset[subset$signals==3,2],alternative="l")
wilcox.test(subset[subset$signals==3,2],subset[subset$signals==4,2],alternative="l")
data2<-data[data$V5>0,]
bin<-c(rep(1,41294),rep(2,41294),rep(3,41294),rep(4,41297))
data<-data[order(data$V6),c(4:6)]
colnames(data)<-c("H3K27ac","H3K27me3","NTR")
data[data$H3K27ac<0,1]=0
data[data$H3K27me3<0,2]=0
logData<-log2(data[,c(1:2)]+0.0001)
data[,1]=rescale(data[,1],newrange = c(0,1))
data[,2]=rescale(data[,2],newrange = c(0,1))
pheatmap(as.matrix(data[,c(1,2)]),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
pheatmap(as.matrix(logData$H3K27ac),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,color = colorRampPalette(c("blue","yellow"))(5),breaks = seq(-4,4,length.out = 6))
##Final figures
data<-read.delim(file="nucleosome.H3K27ac.H3K27me3.bdgcmp.NTR.tsv",header=F)
data<-data[data$V7<0.05 & data$V8>0.8,c(4,5,6)]
data<-data[order(data$V6,decreasing = F),]
bin<-c(rep(1,10997),rep(2,10997),rep(3,10997),rep(4,10999))
boxplot(data$V4~bin,outline=F,ylab="H3K27ac fold enrichment",xlab="Nucleosome turnover rate(low->high)")
boxplot(data$V5~bin,outline=F,ylab="H3K27me3 fold enrichment",xlab="Nucleosome turnover rate(low->high)")
subset<-as.data.frame(cbind(bin,data$V4))
wilcox.test(subset[subset$bin==1,2],subset[subset$bin==2,2],alternative="l")
wilcox.test(subset[subset$bin==2,2],subset[subset$bin==3,2],alternative="l")
wilcox.test(subset[subset$bin==3,2],subset[subset$bin==4,2],alternative="l")
subset<-as.data.frame(cbind(bin,data$V5))
wilcox.test(subset[subset$bin==1,2],subset[subset$bin==2,2],alternative="g")
wilcox.test(subset[subset$bin==2,2],subset[subset$bin==3,2],alternative="g")
wilcox.test(subset[subset$bin==3,2],subset[subset$bin==4,2],alternative="g")

###Revisoin for Molecular Cell
setwd("D:/百度云同步盘/my PC D/projects/help/nucleosome/NewRst/Figures/Fig5/EEDpeaks")
data<-read.delim(file="EED.peaks.repMeansignals.wt.ko.NTR.bed3+",header=F)
data<-as.data.frame(cbind(data$V4,data$V6-data$V5))
colnames(data)<-c("signal","NTRdiff")
data<-data[order(data$signal),]
data$signal<-log2(data$signal)
dcols<-densCols(data$signal, colramp = colorRampPalette(blues9[-(1:3)]))
plot(data$signal,data$NTRdiff,pch=20,col=dcols,xlab="log2(EED signals)",ylab="HTR change(CKO-Het)")
text(12,0,"R-square=0.04")
data<-read.delim(file="EED.peaks.wtH3K27ac.woH3K27me3.signal.wt.ko.NTR.tsv",header=F)
data<-data[order(data$V4),]
bin<-c(rep(1,246),rep(2,246),rep(3,246),rep(4,247))
subset<-as.data.frame(cbind(bin,data$V5))
subset$bin<-as.factor(subset$bin)
ggplot(subset,aes(x=bin, y=V2, fill=bin))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-0.5,1.5)+labs(y="Het HTR",x="EED signals(From low to high)",title="H3K27ac associated EED peaks")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA),plot.title = element_text(hjust = .5))
bin<-c(rep(1,328),rep(2,328),rep(3,328))
setwd("../H3K27ac")
data<-read.delim(file="peaks.ko.wt.signalFC.gt1.5.EEDsignals.H3K27acko_wt.NTRdiff.tsv",header=F)
dcols<-densCols(log2(data$V4), colramp = colorRampPalette(blues9[-(1:3)]))
plot(log2(data$V4),data$V7,pch=20,col=dcols,xlab="log2(EED signals)",ylab="HTR change(CKO-Het)")
data<-read.delim(file="peaks.ko.wt.signalFC.gt1.5.EEDsignals.H3K27acko_wt.NTRdiff.intersectEEDpeak.tsv",header=F)
boxplot(data[data$V9<0,7],data[data$V9>0,7],names=c("EED peak-","EED peak+"),ylab="HTR change(CKO-Het)",xlab="H3K27ac peaks in Figure 5C")
text(1.5,1,"p-value = 3.938e-08")
data<-read.delim(file="peaks.signalFC.gt1.5.EEDsignal.wt.koNTR.tsv",header=F)
data<-data[order(data$V4),]
plot(log2(data$V4),data$V5,pch=20,col=dcols,xlab="log2(EED signals)",ylab="Het HTR",ylim=c(-2,2),xlim=c(-8,6))
bin<-c(rep(1,length(data[data$V4<0,5])),rep(2,length(data[data$V4>=0 & data$V4<1,5])),rep(3,length(data[data$V4>=1 & data$V4<2,5])),rep(4,length(data[data$V4>=2,5])))
subset<-as.data.frame(cbind(bin,data$V5))
subset$bin<-as.factor(subset$bin)
ggplot(subset,aes(x=bin, y=V2, fill=bin))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-2,2)+labs(y="Het HTR",x="EED signals(From low to high)",title="Upregulated H3K27ac peaks in CKO")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA),plot.title = element_text(hjust = .5))
setwd("../H3K27me3")
data<-read.delim(file="peaks.ko.wt.signalFC.lt0.5.EEDsignals.H3K27me3ko_wt.NTRdiff.intersectEEDpeak.tsv",header=F)
boxplot(data[data$V9<0,7],data[data$V9>0,7],names=c("EED peak-","EED peak+"),ylab="HTR change(CKO-Het)",xlab="H3K27me peaks in Figure 5C")
data<-read.delim(file="peaks.signalFC.lt0.5.EEDsignal.wt.koNTR.tsv",header=F)
bin<-c(rep(1,length(data[data$V4<0,5])),rep(2,length(data[data$V4>=0 & data$V4<0.2,5])),rep(3,length(data[data$V4>=0.2 & data$V4<0.5,5])),rep(4,length(data[data$V4>=0.5,5])))
subset<-as.data.frame(cbind(bin,data$V5))
subset$bin<-as.factor(subset$bin)
ggplot(subset,aes(x=bin, y=V2, fill=bin))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-2,2)+labs(y="Het HTR",x="EED signals(From low to high)",title="Upregulated H3K27me3 peaks in CKO")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA),plot.title = element_text(hjust = .5))
bin<-c(rep(1,10244),rep(2,10244),rep(3,10244),rep(4,10244))

setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/Figures/Remodelers")
data<-read.delim(file="swi-snf.fpkm.tsv",header=T,row.names = 1)
data<-read.delim(file="remodeler.fpkm.tsv",header=T,row.names = 1)
logData<-log2(data+1)
logData<-logData[order(logData$EED.Het_1),]
marker<-read.delim(file="markerGene.fpkm.tsv",header = T,row.names = 1)
logMarker<-log2(marker+1)
logData<-rbind(logMarker,logData)
pheatmap(logData,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),breaks = seq(0,8,length.out = 51),cluster_rows = F,cluster_cols = F)
EED<-logData[,c(1:4)]
pheatmap(EED,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),breaks = seq(-2,2,length.out = 51),cluster_rows = F,cluster_cols = F,scale = "row")
EED<-cbind(rowMeans(data[,c(1,2)]),rowMeans(data[,c(3,4)]))
colnames(EED)<-c("Het","EEDCKO")
pheatmap(log2(EED+1),color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_rows = F,cluster_cols = F,scale = "row")

pheatmap(logData,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_rows = T,cluster_cols = F,scale = "row")
logData<-logData[-3,]

#Date 2018-2-7 Brg1 ChIP-seq data 
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/Brg1_inhouse/result")
H3k27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Brg1Signals.diff.HTRdiff.bed3+",header=F)
ggplot(H3k27ac,aes(x=V6, y=V7))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-2,2))
upHTR<-read.delim(file="HTRdown.1kb.down5000.Brg1Signals.diff.HTRdiff.bed3+",header=F)
ggplot(upHTR,aes(x=V6, y=V7))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,-0.5))+scale_x_continuous(limits = c(-1,1))
H3K27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Brg1Signals.diff.H3K27acFC.bed3+",header = F)
ggplot(H3K27ac,aes(x=V4, y=log2(V5)))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="log2(Brg1 signals(CKO-WT))",y="H3K27ac signals(CKO/WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-2,2))
H3K27me3<-read.delim(file="H3K27me3.downRegulated.peaks.Brg1Signals.diff.H3K27me3FC.bed3+",header = F)
ggplot(H3K27me3,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="H3K27me3 signals(CKO/WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,0.5))+scale_x_continuous(limits = c(-2,2))
#Data 2018-2-25 Chd4 Chip-seq data
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/CHD4_inhouse/result")
H3k27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Chd4Signals.diff.HTRdiff.bed3+",header=F)
ggplot(H3k27ac,aes(x=V6, y=V7))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Chd4 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-2,2))+geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

H3K27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Chd4Signals.diff.H3K27acFC.bed3+",header = F)
ggplot(H3K27ac,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Chd4 signals(CKO-WT)",y="H3K27ac signals(CKO/WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(1.5,100))+scale_x_continuous(limits = c(-2,2))+geom_vline(xintercept = 0)
H3K27me3<-read.delim(file="H3K27me3.downRegulated.peaks.Chd4Signals.diff.H3K27me3FC.bed3+",header = F)
ggplot(H3K27me3,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Chd4 signals(CKO-WT)",y="H3K27me3 signals(CKO/WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,0.5))+scale_x_continuous(limits = c(-2,2))+geom_vline(xintercept = 0)
#Date 2018-4-3 Brg1 ChIP-seq data 
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/Brg1_inhouse_080329/lib25/result")
H3k27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Brg1Signals.diff.HTRdiff.bed3+",header=F)
ggplot(H3k27ac,aes(x=V6, y=V7))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-2,2))
H3K27me3<-read.delim(file="H3K27me3.downRegulated.peaks.Brg1Signals.diff.HTRdiff.bed3+",header=F)
ggplot(H3K27me3,aes(x=V6, y=V7))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="NTR(CKO-WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(-20,10))

H3K27ac<-read.delim(file="H3K27ac.upRegulated.peaks.Brg1Signals.diff.H3K27acFC.bed3+",header = F)
ggplot(H3K27ac,aes(x=V4, y=log2(V5)))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="log2(H3K27ac signals(CKO/WT))")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
H3K27me3<-read.delim(file="H3K27me3.downRegulated.peaks.Brg1Signals.diff.H3K27me3FC.bed3+",header = F)
ggplot(H3K27me3,aes(x=V4, y=V5))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+labs(x="Brg1 signals(CKO-WT)",y="H3K27me3 signals(CKO/WT)")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,0.5))+scale_x_continuous(limits = c(-20,10))

#Final figure 6
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/H2BGFP-Brg1KD-week0-old180403/mergedBam")
data<-read.delim(file="Brg1KD.scramble.H2BGFP.tsv",header=F)
colnames(data)<-c("Brg1KD","scramble","region")
data$Brg1KD<-log2(data$Brg1KD)
data$scramble<-log2(data$scramble)
levels(data$region)<-c("H3K27me3","H3K27ac","EED")
ggplot(data,aes(x=Brg1KD,y=scramble))+geom_point(size=2,shape=21,aes(col=region))+ylim(-6,10)+xlim(-6,10)+geom_abline(slope = 1,intercept = 0)+coord_equal(ratio = 1)
data2<-melt(data)
ggplot(data2,aes(x=region,y=value,col=variable))+geom_boxplot(outlier.shape = NA)+labs(x="",y="log2(H2BGFP signals)")+ylim(-5,2.5)
EED<-subset(data,region=="EED")
ggplot(EED,aes(x=Brg1KD, y=scramble))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+ylim(-6,10)+xlim(-6,10)+geom_abline(slope = 1,intercept = 0)+coord_equal(ratio = 1)+labs(x="log2(Brg1kd)",y="log2(scramble)",title="EED")
H3K27ac<-subset(data,region=="H3K27ac")
ggplot(H3K27ac,aes(x=Brg1KD, y=scramble))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+ylim(-6,10)+xlim(-6,10)+geom_abline(slope = 1,intercept = 0)+coord_equal(ratio = 1)+labs(x="log2(Brg1kd)",y="log2(scramble)",title="H3K27ac")
H3K27me3<-subset(data,region=="H3K27me3")
ggplot(H3K27me3,aes(x=Brg1KD, y=scramble))+stat_bin_2d(bins=200)+scale_fill_gradient(low="lightblue",high="red")+ylim(-6,10)+xlim(-6,10)+geom_abline(slope = 1,intercept = 0)+coord_equal(ratio = 1)+labs(x="log2(Brg1kd)",y="log2(scramble)",title="H3K27me3")
data3<-cbind(data[,1]-data[,2],as.character(data[,3]))
data3<-as.data.frame(data3,stringsAsFactors =F)
data3$V2<-as.factor(data3$V2)
data3$V1<-as.numeric(data3$V1)
data3<-na.omit(data3)
data3<-data3[is.finite(data3$V1),]
library("plyr")
mu<-ddply(data3,.(V2),summarise, grp.mean=mean(V1))
ggplot(data3,aes(x=V1,color=V2))+geom_density()+scale_x_continuous(limits = c(-5,5))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="log2(H2BGFP signals(BRG1KD/WT))",y="Estimated density of peak count")
ggplot(data3,aes(x=V1,..count..,color=V2))+geom_density()+scale_x_continuous(limits = c(-5,5))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="log2(H2BGFP signals(BRG1KD/WT))",y="Peak count")
#barplot
x<-matrix(c(2952,190,6022,7857,440,19035,31141,1031,20419),nrow=3)
rownames(x)<-c("decrease","unchange","increase")
colnames(x)<-c("EED","H3K27ac","H3K27me3")
x[,1]=x[,1]/sum(x[,1])
x[,2]=x[,2]/sum(x[,2])
x[,3]=x[,3]/sum(x[,3])
barplot(x,ylab="Fraction",col=c(gray(0.1),gray(0.4),gray(0.9)),ylim=c(0,1.3))
legend("topright",c("increase","unchange","decrease"),fill=c(gray(0.9),gray(0.4),gray(0.1)))
x<-round(c(x[1,1],x[2,1],x[1,2],x[2,2]),digits = 4)*100
percent=""
percent[1]=paste(x[1],"%",sep="")
percent[2]=paste(x[2],"%",sep="")
percent[3]=paste(x[3],"%",sep="")
percent[4]=paste(x[4],"%",sep="")
text(c(0.8,0.8,1.8,1.8),c(0.3,0.8,0.3,0.8),percent)
fisher.test(x) #p-value <2.2e-16


##Add NTR change for EED peaks
data<-read.delim(file="D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/Figures/Fig5/EED/EEDpeaks.ko-wt.NTR.diff.EEDsignals.tsv",header=F)
ggplot(data,aes(x=V7, y=V6))+stat_bin_2d(bins=200)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,20),na.value = "red")+labs(x="EED signals",y="NTR(CKO-WT)")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(-2,2))+scale_x_continuous(limits = c(0,10))

##Add cross-analysis for EED and Brg1 peaks 
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/Brg1-ChIP-lib1825-Work/Brg1-ChIP/result/EED_classify")
wt<-read.delim(file="EED.wtBrg1Binding.wt.ko.NTR.tsv",header = F)
wo<-read.delim(file="EED.woBrg1Binding.wt.ko.NTR.tsv",header=F)
subset<-rbind(cbind(wt$V5-wt$V4,"wtBrg1"),cbind(wo$V5-wo$V4,"woBrg1"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
ggplot(subset,aes(x=V2, y=V1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA,color=c("blue","red"))+
  ylim(-1,0.8)+labs(y="HTR change(CKO-Het)",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA),plot.title = element_text(hjust = .5))

library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,color=V2))+geom_density()+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KO-WT)")
pheatmap(subset[subset$V2=="wtBrg1",],cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(10),breaks=seq(-1,1,length.out = 11))
pheatmap(wo$V5-wo$V4,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(10),breaks=seq(-1,1,length.out = 11))

data<-read.delim(file="EEDpeaks.Brg1Signals.sorted.wt.ko.NTR.tsv",header=F)
pheatmap(data$V5-data$V4,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(10),breaks=seq(-1,0.5,length.out = 11))
pheatmap(cbind(data$V5,data$V4),cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(10),breaks=seq(-0.5,0.5,length.out = 11))

agg<-read.delim(file="EED.woBrg1Binding.NC.agg.txt",header=F)
agg<-read.delim(file="EED.wtBrg1Binding.NC.agg.txt",header=F)
par(mfrow=c(2,1))
data<-as.matrix(agg[,c(2:5)])
#Norm by whole region mean(Final version)
rst<-data/t(replicate(10000,colMeans(data)))
plot(agg[,1],(rst[,1]+rst[,2])/2,type="l",col="red",ylim=c(0.6,1.2),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
lines(agg[,1],(rst[,3]+rst[,4])/2,col="blue")
legend("topright",legend = c("WT","EEDko"),col=c("red","blue"),lty=1)

agg<-read.delim(file="EED.woBrg1Binding.Brg1KD.agg.txt",header=F)
agg<-read.delim(file="EED.wtBrg1Binding.Brg1KD.agg.txt",header=F)
par(mfrow=c(2,1))
data<-as.matrix(agg[,c(2:3)])
#Norm by whole region mean(Final version)
rst<-data/t(replicate(10000,colMeans(data)))
plot(agg[,1],rst[,1],type="l",col="red",ylim=c(0.6,1.2),xlab="Peak center relative position",ylab="Normalized nucleosome occupancy")
lines(agg[,1],rst[,2],col="blue")
legend("topright",legend = c("WT","Brg1KD"),col=c("red","blue"),lty=1)

#20180914 new data processing
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/H2BGFP-Brg1KD-week0-20180914/mergedBam/NTR")
files=list.files("Brg1KD.NTR.tsv$")
for(i in 1:length(files)){
  controlFile=sub(".Brg1KD.NTR.tsv",".Control.NTR.tsv",files[i])
  prefix=sub(".Brg1KD.NTR.tsv","",files[i])
  brg1kd=read.delim(file=files[i],header=F)
  control=read.delim(file=controlFile,header = F)
  subset<-cbind(brg1kd$V6,control$V6)
  colnames(subset)<-c("Brg1KD","Control")
  subset<-melt(subset,value.name = "HTR")
  pdf(file=paste(prefix,"NTR.pdf",sep=""),width = 4,height = 8)
  ggplot(subset,aes(x=Var2, y=HTR, fill=Var2))+geom_boxplot()+labs(y="HTR",x="")+
    stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
  #boxplot(brg1kd$V6,control$V6,names = c("Brg1KD","Control"),ylab="HTR",outline = F,main=prefix,col="gray")
  dev.off()
}

#NTR change in Brg1kd
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/ChIP-seq/others/H2BGFP-Brg1KD-week0-20180914/result/NTR")
h3k27ac<-read.delim(file="H3K27ac.upRegulated.peaks.kd.wt.NTR.tsv",header=F)
h3k27me3<-read.delim(file="H3K27me3.downRegulated.peaks.kd.wt.NTR.tsv",header=F)
subset<-rbind(cbind(h3k27ac$V11-h3k27ac$V12,"H3K27ac"),cbind(h3k27me3$V11-h3k27me3$V12,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,color=V2))+geom_density()+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KD-WT)")
ks.test(subset[subset$V2=="H3K27ac",1],subset[subset$V2=="H3K27me3",1])#p-value < 2.2e-16

plot(ecdf(subset[subset$V2=="H3K27ac",1]),xlim=c(-2,2),verticals = T, do.points=F,xlab="NTR(KD-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(subset[subset$V2=="H3K27me3",1]),col="red",verticals = T, do.points=F)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("black","red"), lty=1)
#NTR in EED-Brg1 co-bind regions
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/revision_201809")
EED<-read.delim(file="EED.unique.NTR.tsv",header=F)
common<-read.delim(file="EED.Brg1.commonReion.NTR.tsv",header=F)
subset<-rbind(cbind(common$V10,"wtBrg1"),cbind(EED$V4,"woBrg1"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
ggplot(subset,aes(x=V2, y=V1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA,color=c("blue","red"))+
  ylim(-0.5,1.5)+labs(y="HTR",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA),plot.title = element_text(hjust = .5))
#Circulation research revision
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/Fig3D")
random_ac<-read.delim(file="H3K27ac.upRegulated.peaks.shuffle.wt.ko.NTRdiff.tsv",header=F)
random_me3<-read.delim(file="H3K27me3.downRegulated.peaks.shuffle.wt.ko.NTRdiff.tsv",header=F)
ac<-read.delim(file="../../Figures/Fig5/H3K27ac/peaks.ko.wt.signalFC.gt1.5.EEDsignals.H3K27acko_wt.NTRdiff.tsv",header = F)
me3<-read.delim(file="../../Figures/Fig5/H3K27me3/peaks.ko.wt.signalFC.lt0.5.EEDsignals.H3K27me3ko_wt.NTRdiff.tsv",header=F)
subset<-rbind(cbind(ac$V7,"H3K27ac"),cbind(me3$V7,"H3K27me3"),cbind(random_ac$V6,"random_H3K27ac"),cbind(random_me3$V6,"random_H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+scale_fill_manual(values=c("#D4A6F7", "#ADC879", "#F1A6A4", "#82D1D3"))+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KO-WT)")
plot(ecdf(me3$V7),xlim=c(-1,1),verticals = T, col="#D4A6F7", do.points=F,xlab="NTR(CKO-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(random_me3$V6),col="#ADC879",verticals = T, do.points=F)
abline(v=0,col="gray",lty=2)
abline(h=0.5,col="gray",lty=2)
legend("topleft",c("H3K27me3","Random_H3K27me3"),col=c("#D4A6F7","#ADC879"), lty=1)

subset<-rbind(cbind(random_ac$V6,"random_H3K27ac"),cbind(random_me3$V6,"random_H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+scale_fill_manual(values=c("#D4A6F7", "#ADC879"))+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KO-WT)")
plot(ecdf(random_ac$V6),xlim=c(-1,1),verticals = T, col="#D4A6F7", do.points=F,xlab="NTR(CKO-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(random_me3$V6),col="#ADC879",verticals = T, do.points=F)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("#D4A6F7","#ADC879"), lty=1)


allac<-read.delim(file="H3K27ac-WT.peaks.wt.ko.NTRdiff.tsv",header=F)
allme3<-read.delim(file="H3K27me3-WT.peaks.wt.ko.NTRdiff.tsv",header=F)
subset<-rbind(cbind(allac$V6,"H3K27ac"),cbind(allme3$V6,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(KO-WT)")

plot(ecdf(allac$V6),xlim=c(-1,1),verticals = T, do.points=F,xlab="NTR(CKO-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(allme3$V6),col="red",verticals = T, do.points=F)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("black","red"), lty=1)

H3K27ac<-read.delim(file="H3K27ac.up.ko-wtNCdiff.tsv",header=F)
H3K27me3<-read.delim(file="H3K27me3.down.kp-wtNCdiff.tsv",header=F)
subset<-rbind(cbind(H3K27ac$V6,"H3K27ac"),cbind(H3K27me3$V6,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
ggplot(subset,aes(x=V2, y=V1))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-0.5,0.5)+labs(y="Nucleosome occupancy change (CKO-Het)",x="")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(H3K27ac$V6,H3K27me3$V6)#0.0003

agg<-read.delim(file="H3K27ac.up.EEDHet.ko.withCI.agg.tsv",header=F)
agg<-read.delim(file="H3K27me3.down.EEDHet.ko.withCI.agg.tsv",header=F)
data<-as.matrix(agg[,c(2:7)])
#rst<-data/t(replicate(10000,colMeans(data)))
rst<-data/t(replicate(10000,c(rep(colMeans(data)[1],3),rep(colMeans(data)[4],3))))
data<-rbind(cbind(class=replicate(10000,"Het"),pos=agg[,1],Mean=rst[,1],Low=rst[,2],High=rst[,3]),cbind(class=replicate(10000,"CKO"),pos=agg[,1],Mean=rst[,4],Low=rst[,5],High=rst[,6]))
data<-as.data.frame(data,stringsAsFactors = F,row.names = NULL)
data$class<-as.factor(data$class)
data$pos<-as.numeric(data$pos)
data$Mean<-as.numeric(data$Mean)
data$Low<-as.numeric(data$Low)
data$High<-as.numeric(data$High)
ggplot(data, aes(pos, Mean))+geom_ribbon(aes(ymin = Low, ymax = High, fill = class), alpha = .25) + geom_line(aes(colour = class)) +
  labs(x = "Distance to peak center(bp)", y = "Nucleosome density")+ylim(0.9,1.1)+scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
subset<-data[data$pos>=-1500 & data$pos<=1500,]
subset<-data[data$pos>=-2500 & data$pos<=2500,]
subset<-data[data$pos>=-500 & data$pos<=500,]
ggplot(subset, aes(pos, Mean))+geom_ribbon(aes(ymin = Low, ymax = High, fill = class), alpha = .25) + geom_line(aes(colour = class)) +
  labs(x = "Distance to peak center(bp)", y = "Nucleosome density")+ylim(0.98,1.08)+scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/TF_control")
data<-read.delim(file="Klf4.Nanog.Oct4.Sox2.txt",header=F)
data<-data[101:401,]
data_f<-read.delim(file="Klf4.Nanog.Oct4.Sox2.rmHeartTF.txt",header=F)
data_f<-data_f[101:401,]
par(mfrow=c(1,2))
myColors=rainbow(10)[1:3]
plot(seq(-3000,3000,length.out = 301),data[,2],xlab="Distance to peak center(bp)",ylab="HTR",type="l",ylim=c(0,0.6),col=myColors[1],main="All peaks")
for(i in 3:4){
  lines(seq(-3000,3000,length.out = 301),data[,i],col=myColors[i-1])
}
legend("topright",c("Nanog","Oct4","Sox2"),col=myColors,lty=1)
plot(seq(-3000,3000,length.out = 301),data_f[,2],xlab="Distance to peak center(bp)",ylab="HTR",type="l",ylim=c(0,0.6),col=myColors[1],main="Peaks (rm Heart TF)")
for(i in 3:4){
  lines(seq(-3000,3000,length.out = 301),data_f[,i],col=myColors[i-1])
}
legend("topright",c("Nanog","Oct4","Sox2"),col=myColors,lty=1)
enhancer<-read.delim(file="enhancers.TFs.status.NTR.tsv",header=F)
data<-matrix(nrow=2,ncol=3)
data[,1]<-c("enhancers","enhancers")
data[,2]<-c("wt","wo")
data[,3]<-c(0.2191,0.7809)#wt:6269,wo:22337
data<-as.data.frame(data,stringsAsFactors=F)
colnames(data)=c("class","TFs","Fraction")
data$TFs<-factor(data$TFs,levels=c("wt","wo"))
data$Fraction<-round(as.numeric(data$Fraction),digits = 2)
ggplot(data,aes(x=class,y=Fraction,fill=TFs))+geom_bar(stat="identity")+labs(y="Fraction")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
data<-as.data.frame(enhancer[,c(4,5)])
ggplot(data,aes(x=V4, y=V5, fill=V4))+geom_boxplot(width=0.4,position=position_dodge(0.9))+ylim(-0.5,1.0)+labs(y="NTR")+stat_summary(fun.y=mean,position = position_dodge(width = .9),geom="point", size=3,color="white")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(enhancer[enhancer$V4=="wt",5],enhancer[enhancer$V4=="wo",5])#p=0.2859
enhancer<-read.delim(file="enhancers.TFnum.NTR.final.tsv",header=F)
enhancer<-enhancer[enhancer$V3-enhancer$V2>20,]
enhancer$V4<-enhancer$V4-1
subset<-as.data.frame(enhancer[,c(4,5)])
colnames(subset)<-c("status","NTR")
subset$status<-factor(subset$status)
ggplot(subset,aes(x=status, y=NTR, fill=status))+geom_boxplot(lwd=1,fill="gray")+labs(x="TFs binding number",y="NTR")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  stat_summary(fun.y=mean, geom="point", size=2,color="white")
wilcox.test(subset[subset$status=="1",2],subset[subset$status=="2",2],alternative = "l")#0.007
wilcox.test(subset[subset$status=="2",2],subset[subset$status=="3",2],alternative = "l")#0.99

data<-read.delim(file="enhancers.TFnum.allTFsignal.NTR.tsv",header=F)
data<-data[data$V3-data$V2>20,]
data<-na.omit(data)
data[data$V4<0,4]=0
data[data$V5<0,5]=0
data[data$V6<0,6]=0
data<-data[order(rowMeans(data[,c(4:6)]),decreasing = F),]
bin<-c(rep(1,8036),rep(2,910),rep(3,574))#[0,10],(10,20),[20,)]
subset<-as.data.frame(cbind(bin,data$V8))
colnames(subset)<-c("Signals","NTR")
subset$Signals<-factor(subset$Signals)
ggplot(subset,aes(x=Signals, y=NTR, fill=Signals))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.2,lwd=1,fill=NA)+
  ylim(-1,1.6)+labs(y="NTR",x="TF signals")+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
wilcox.test(subset[subset$Signals=="1",2],subset[subset$Signals=="2",2],alternative = "l")#<2.2e-16; 2.822e-09

setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/diffIntervals")
library("RColorBrewer")
data<-read.delim(file="mm10.3kbIntervals.Anno.summary.tsv",header=F)
random<-read.delim(file="random3k.Anno.summary.tsv",header=F)
fisher.test(matrix(c(data[1,2],10000-data[1,2],random[1,2],10000-random[1,2]),byrow = F,nrow = 2))
data<-cbind(data$V2,random$V2)
colnames(data)<-c("Top10k","random10k")
rownames(data)<-random$V1
data[,1]<-data[,1]/sum(data[,1])
data[,2]<-data[,2]/sum(data[,2])
barplot(data,horiz=T,col=brewer.pal(n = 12, name = "Paired")[c(2,1,3,6,7,10,12)],legend=rownames(data),xlab="Percentage")
text(c(0.05,0.05),c(1,2),c(round(data[1,1],digits = 4),round(data[1,2],digits = 4)))
data<-read.delim(file="mm10.2kbIntervals.Anno.summary.tsv",header=F)
random<-read.delim(file="random2k.Anno.summary.tsv",header=F)

data<-read.delim(file="mm10.500bpIntervals.Anno.summary.tsv",header=F)
random<-read.delim(file="random500bp.Anno.summary.tsv",header=F)

#New GOF data (20190227)
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/gainOfFunction/EEDGOF_P1/result")
H3K27ac<-read.delim(file="H3K27ac.up.HTR.control.GOF.GOF-Control.tsv",header=F)
H3K27me3<-read.delim(file="H3K27me3.down.HTR.control.GOF.GOF-Control.tsv",header = F)
subset<-rbind(cbind(H3K27ac$V6,"H3K27ac"),cbind(H3K27me3$V6,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
ggplot(subset,aes(x=V2,y=V1,fill=V2))+geom_violin(alpha=0.5,trim=F)+geom_boxplot(width=0.4, lwd=1,fill=NA,color="black")+scale_y_continuous(limits = c(-1.5,1.5))+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="",y="HTR(GOF-WT)")
x<-matrix(c(nrow(H3K27ac[H3K27ac$V6<0,]),nrow(H3K27me3[H3K27me3$V6<0,]),nrow(H3K27ac[H3K27ac$V6>=0,]),nrow(H3K27me3[H3K27me3$V6>=0,])),nrow=2,byrow=T)
rownames(x)<-c("decrease","increase")
colnames(x)<-c("H3K27ac","H3K27me3")
barplot(x/matrix(c(colSums(x),colSums(x)),nrow=2,byrow=T),ylab="Peak number",col=c(gray(0.4),gray(0.9)))
legend("topright",c("increase","decrease"),fill=c(gray(0.9),gray(0.4)))
#H3K27ac  H3K27me3
#increase 0.3652225 0.4508249
#decrease 0.6347775 0.5491751
plot(ecdf(H3K27ac$V6),xlim=c(-1,1),verticals = T, do.points=F,xlab="NTR(GOF-WT)",ylab="Cumulative frequency",main="")
lines(ecdf(H3K27me3$V5),col="red",verticals = T, do.points=F)
legend("topleft",c("H3K27ac","H3K27me3"),col=c("black","red"), lty=1)

library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(GOF-WT)")
files<-list.files(pattern="H3K27ac.*control.GOF.HTR.tsv$")
for(i in 1:length(files)){
  me3File=sub("^H3K27ac","H3K27me3",files[i])
  prefix=paste("H3K27me3.",sub("tsv","pdf",files[i]))
  H3K27ac<-read.delim(file=files[i],header=F)
  H3K27me3<-read.delim(file=me3File,header = F)
  subset<-rbind(cbind(H3K27ac$V6,"H3K27ac"),cbind(H3K27me3$V6,"H3K27me3"))
  subset<-as.data.frame(subset,stringsAsFactors =F)
  subset$V1<-as.numeric(subset$V1)
  subset$V2<-as.factor(subset$V2)
  mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
  pvalue<-ks.test(H3K27ac$V6,H3K27me3$V6)$p.value
  if(pvalue == 0){pvalue="p < 2.2e-16"}else{pvalue=paste("p=",pvalue,sep="")}
  pdf(file=prefix)
  ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="HTR(GOF-WT)")+annotate("text",x=0.5,y=2,label=pvalue)
  dev.off()
}
eed<-read.delim(file="EED.peaks.control.GOF.HTR.tsv",header = F)
boxplot(eed$V4,eed$V5,outline = F,names = c("Control","EEDGOF"),ylab="HTR")
#New replicates 
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/gainOfFunction/replicates_supp/result")
data<-read.delim(file="D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/Figures/Fig1/H3.H2BGFP.1kb.sum.tsv",header=F)
data<-as.matrix(data[,2])
new<-read.delim(file="H3.Rep2.1kb.new.simple.tsv",header = F)
data<-na.omit(cbind(new,data))
cor(data[,1],data[,2],method="spearman")# -0.10
cor(data[,1],data[,2])# 0.05
ggplot(data,aes(x=V1, y=data))+stat_bin_2d(bins=100)+scale_fill_gradientn(colours = c("blue","yellow","red"),limits=c(0,2000),na.value="red")+labs(x="H3 signal",y="week0 H2BGFP signal")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+
  scale_y_continuous(limits=c(0,3))+scale_x_continuous(limits = c(0,3))+theme(legend.position=c(0.9,0.9), legend.justification=c(0.3,1))

H3K27ac<-read.delim(file="H3K27ac.up.ko.wt.NTR.diff.tsv",header=F)
H3K27me3<-read.delim(file="H3K27me3.down.ko.wt.NTR.diff.tsv",header = F)
subset<-rbind(cbind(H3K27ac$V6,"H3K27ac"),cbind(H3K27me3$V6,"H3K27me3"))
subset<-as.data.frame(subset,stringsAsFactors =F)
subset$V1<-as.numeric(subset$V1)
subset$V2<-as.factor(subset$V2)
library("plyr")
mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="NTR(CKO-WT)")
files<-list.files(pattern="H3K27ac.*ko.wt.NTR.diff.tsv$")
for(i in 1:3){
  me3File=sub("^H3K27ac","H3K27me3",files[i])
  prefix=paste("H3K27me3.",sub("tsv","pdf",files[i]))
  H3K27ac<-read.delim(file=files[i],header=F)
  H3K27me3<-read.delim(file=me3File,header = F)
  subset<-rbind(cbind(H3K27ac$V6,"H3K27ac"),cbind(H3K27me3$V6,"H3K27me3"))
  subset<-as.data.frame(subset,stringsAsFactors =F)
  subset$V1<-as.numeric(subset$V1)
  subset$V2<-as.factor(subset$V2)
  mu<-ddply(subset,.(V2),summarise, grp.mean=mean(V1))
  pvalue<-ks.test(H3K27ac$V6,H3K27me3$V6)$p.value
  if(pvalue == 0){pvalue="p < 2.2e-16"}else{pvalue=paste("p=",pvalue,sep="")}
  pdf(file=prefix)
  ggplot(subset,aes(x=V1,fill=V2))+geom_density(alpha=0.5)+scale_x_continuous(limits = c(-1,1))+geom_vline(data=mu,aes(xintercept=grp.mean, color=V2),linetype="dashed")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))+labs(x="HTR(CKO-WT)")+annotate("text",x=0.5,y=2,label=pvalue)
  dev.off()
}
eed<-read.delim(file="EED.peaks.ko.wt.NTR.diff.tsv",header = F)
boxplot(eed$V5,eed$V4,outline = F,names = c("Control","EEDko"),ylab="HTR")

setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/CirRes_revision201812/gainOfFunction/MNase_20190309/result/")
agg<-read.delim(file="H3K27ac.up.Control.EEDGOF.withCI.agg.tsv",header=F)
agg<-read.delim(file="H3K27me3.down.Control.EEDGOF.withCI.agg.tsv",header=F)
data<-as.matrix(agg[,c(2:7)])
rst<-data/t(replicate(10000,c(rep(colMeans(data)[1],3),rep(colMeans(data)[4],3))))
data<-rbind(cbind(class=replicate(10000,"Het"),pos=agg[,1],Mean=rst[,1],Low=rst[,2],High=rst[,3]),cbind(class=replicate(10000,"CKO"),pos=agg[,1],Mean=rst[,4],Low=rst[,5],High=rst[,6]))
data<-as.data.frame(data,stringsAsFactors = F,row.names = NULL)
data$class<-as.factor(data$class)
data$pos<-as.numeric(data$pos)
data$Mean<-as.numeric(data$Mean)
data$Low<-as.numeric(data$Low)
data$High<-as.numeric(data$High)
ggplot(data, aes(pos, Mean))+geom_ribbon(aes(ymin = Low, ymax = High, fill = class), alpha = .25) + geom_line(aes(colour = class)) +
  labs(x = "Distance to peak center(bp)", y = "Nucleosome density")+ylim(0.9,1.1)+scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
setwd("D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/MNase-seq/MNaseSeq-Brg1KD-week0-20180914/result")
agg<-read.delim(file="H3K27ac.up.Scramble.Brg1KD.withCI.agg.tsv",header=F)
agg<-read.delim(file="H3K27me3.down.Scramble.Brg1KD.withCI.agg.tsv",header=F)
data<-as.matrix(agg[,c(2:7)])
rst<-data/t(replicate(10000,c(rep(colMeans(data)[1],3),rep(colMeans(data)[4],3))))
data<-rbind(cbind(class=replicate(10000,"Scramble"),pos=agg[,1],Mean=rst[,1],Low=rst[,2],High=rst[,3]),cbind(class=replicate(10000,"BRG1KD"),pos=agg[,1],Mean=rst[,4],Low=rst[,5],High=rst[,6]))
data<-as.data.frame(data,stringsAsFactors = F,row.names = NULL)
data$class<-as.factor(data$class)
data$pos<-as.numeric(data$pos)
data$Mean<-as.numeric(data$Mean)
data$Low<-as.numeric(data$Low)
data$High<-as.numeric(data$High)
ggplot(data, aes(pos, Mean))+geom_ribbon(aes(ymin = Low, ymax = High, fill = class), alpha = .25) + geom_line(aes(colour = class)) +
  labs(x = "Distance to peak center(bp)", y = "Nucleosome density")+ylim(0.9,1.1)+scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
subset<-data[data$pos>=-2500 & data$pos<=2500,]
subset<-data[data$pos>=-1500 & data$pos<=1500,]
subset<-data[data$pos>=-500 & data$pos<=500,]
ggplot(subset, aes(pos, Mean))+geom_ribbon(aes(ymin = Low, ymax = High, fill = class), alpha = .25) + geom_line(aes(colour = class)) +
  labs(x = "Distance to peak center(bp)", y = "Nucleosome density")+ylim(0.9,1.1)+scale_color_manual(values=c("blue","red"))+scale_fill_manual(values=c("blue","red"))+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))

#20190416 Revision2
data<-read.delim(file="D:/百度云同步盘/my PC D/projects/nucleosomeTurnOver/NewRst/Figures/Fig5/H3K27me3/peaks.ko.wt.signalFC.lt0.5.wt_ko_diffHTR.tsv",header=F)
wilcox.test(data$V4,data$V5,paired = T)#p-value < 2.2e-16
colnames(data)[4]<-"Het"
colnames(data)[5]<-"CKO"
data<-melt(data[,c(4,5)])
data$variable<-factor(data$variable)
ggplot(data,aes(x=variable, y=value))+geom_violin(trim=F,fill="gray",color=NA)+geom_boxplot(width=0.4, lwd=1,fill=NA)+labs(y="HTR",x="")+
  ylim(-1,1)+stat_summary(fun.y=mean, geom="point", size=2,color="white")+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="white"),panel.border = element_rect(colour="black",fill=NA))
