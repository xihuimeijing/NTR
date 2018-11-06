setwd("F:/Project/OKCAM")
library(VennDiagram)
###################################################################################
#Figure 1
###################################################################################
InterPro=read.table("figure/data/CAM.venn.InterPro",header=F)
InterPro<-InterPro[,1]
GO=read.table("figure/data/CAM.venn.GO",header=F)
GO<-GO[,1]
George=read.table("figure/data/CAM.venn.George",header=F)
George<-George[,1]
#Length_InterPro<-length(InterPro)
#Length_GO<-length(GO)
#Length_George<-length(George)
#Length_InterProGO<-length(intersect(InterPro,GO))
#Length_InterProGeorge<-length(intersect(InterPro,George))
#Length_InterProGOGeorge<-length(intersect(intersect(InterPro,George),GO))
#Length_GOGeorge<-length(intersect(GO,George))
T<-venn.diagram(list(A=InterPro,B=GO,C=George),filename=NULL,
                lwd=1,lty=1,col=c('red','green','blue'),
                alpha=0.3,cex=3,label.col="white",
                fill=c('red','green','blue'),
                cat.col=c('red','green','blue'),
                category.names=c("InterPro","GO","EntrezID Gene Query"),
                cat.cex=2,cat.pos=c(180,-30,30),cat.dist=c(0.03,0.05,0.05),
                rotation.degree=90,)
grid.draw(T)

#Following is OKCAM.v2 specific
InterPro=read.table("figure/Figure1/data/OKCAM.v2.specific.InterPro",header=F)
InterPro<-InterPro[,1]
GO=read.table("figure/Figure1/data/OKCAM.v2.specific.GO",header=F)
GO<-GO[,1]
George=read.table("figure/Figure1/data/OKCAM.v2.specific.George",header=F)
George<-George[,1]
T<-venn.diagram(list(A=InterPro,B=GO,C=George),filename=NULL,
                lwd=1,lty=1,col=c('red','green','blue'),
                alpha=0.3,cex=3,label.col="white",
                fill=c('red','green','blue'),
                cat.col=c('red','green','blue'),
                category.names=c("InterPro","GO","EntrezID Gene Query"),
                cat.cex=2,cat.pos=c(230,-30,30),cat.dist=c(0.03,0.05,0.05),
                rotation.degree=90,)
grid.draw(T)

InterPro=read.table("figure/Figure1/data/OKCAM.v2.specific.InterPro",header=F)
InterPro<-InterPro[,1]
GO=read.table("figure/Figure1/data/OKCAM.v2.specific.GO",header=F)
GO<-GO[,1]
George=read.table("figure/Figure1/data/OKCAM.v2.specific.George",header=F)
George<-George[,1]
T<-venn.diagram(list(A=InterPro,B=GO,C=George),filename=NULL,
                lwd=2,lty=1,col=c('red','black','blue'),
                alpha=0.3,cex=3,
                cat.col=c('red','black','blue'),
                category.names=c("InterPro","GO","EntrezID Gene Query"),
                cat.cex=2,cat.pos=c(230,-30,30),cat.dist=c(0.03,0.05,0.05),
                rotation.degree=90,)
grid.draw(T)

#Following is OKCAM.v1 overlap with OKCAM.v2
OKCAM.v1=read.table("figure/Figure1/data/CAM.OKCAM.v1.entrezID",header=F)
OKCAM.v1<-OKCAM.v1[,1]
OKCAM.v2=read.table("figure/Figure1/data/CAM.OKCAM.v2.entrezID",header=F)
OKCAM.v2<-OKCAM.v2[,1]
#Lenghth_OKCAM.v1<-length(OKCAM.v1)
#Length_OKCAM.v2<-length(OKCAM.v2)
#Length_v1v2<-length(intersect(OKCAM.v1,OKCAM.v2))
#T<-venn.diagram(list(A=OKCAM.v1,B=OKCAM.v2),filename=NULL,
#                lwd=1,col=c("red","blue"),cex=3,alpha=0.3,
#                label.col=c("white","white","white"),
#                fill=c("red","blue"),
#                cat.col=c("red","blue"),
#                cat.pos=c(-30,30),cat.dist=c(0.05,0.05),cat.cex=2,
#                category.names=c("OKCAM v1.0","OKCAM v2.0")
#)
T<-venn.diagram(list(A=OKCAM.v1,B=OKCAM.v2),filename=NULL,
                lwd=2,col=c("red","blue"),cex=3,lty=1,
                cat.pos=c(-30,30),cat.dist=c(0.05,0.05),cat.cex=2,
                rotation.degree=0
)
grid.draw(T)

#FOllowing is OKCAM.v1 specific classcify
dat=read.table("figure/data/OKCAM.v1.specific.classcify",header=F)
ratio=sprintf("%.2f",100*dat[,2]/sum(dat[,2]))
ratio=paste(ratio,"%",sep="")
label=paste(dat[,1],ratio,sep="\n")
pie(dat[,2],col=c("red","green","blue","purple"),border="purple",labels=label,font=2)

dat=read.table("figure/Figure1/data/OKCAM.v1.specific.classcify",header=F)
ratio=sprintf("%.2f",100*dat[,2]/sum(dat[,2]))
ratio=paste(ratio,"%",sep="")
label=paste(dat[,1],ratio,sep="\n")
pie(dat[,2],col=c("red","black","blue"),labels=label,font=2)

OKCAM.add.InterPro=read.table("Figure/data/OKCAM.add.InterPro",header=F)
OKCAM.add.InterPro<-OKCAM.add.InterPro[,1]
OKCAM.add.GO=read.table("Figure/data/OKCAM.add.GO",header=F)
OKCAM.add.GO<-OKCAM.add.GO[,1]
OKCAM.add.George=read.table("Figure/data/OKCAM.add.George",header=F)
OKCAM.add.George<-OKCAM.add.George[,1]
T<-venn.diagram(list(A=OKCAM.add.InterPro,B=OKCAM.add.GO,C=OKCAM.add.George),
                filename=NULL,
                lwd=1,lty=2,col=c("red","green","blue"),
                fill=c('red','green','blue'),
                cat.col=c('red','green','blue'),
                rotation.degree=90)
grid.draw(T)

###################################################################################
#Figure 2
###################################################################################

#Following is alleleRegulation overlap with alleleExpression
regulation=read.table("figure/data/CAM_allRegulate.entrezID",header=F)
regulation<-regulation[,1]
expression=read.table("figure/data/CAM_expression.entrezID",header=F)
expression<-expression[,1]
T<-venn.diagram(list(A=regulation,B=expression),filename=NULL,
                lwd=2,col=c("red","blue"),cex=3,lty=1,
                cat.pos=c(-30,30),cat.dist=c(0.05,0.05),cat.cex=2,
                rotation.degree=0
)
grid.draw(T)

#Following is overlap of each regulation
miRNA=read.table("figure/data/CAM_miRNATarget.entrezID",header=F)
miRNA<-miRNA[,1]
tfbs=read.table("figure/data/CAM_tfbs.entrezID",header=F)
tfbs<-tfbs[,1]
chiapet=read.table("figure/data/CAM_chiapet.entrezID",header=F)
chiapet<-chiapet[,1]
splicing=read.table("figure/data/CAM_splicing.entrezID",header=F)
splicing<-splicing[,1]
other=read.table("figure/data/CAM_otherRegulate.entrezID",header=F)
other<-other[,1]
T<-venn.diagram(list(A=miRNA,B=tfbs,C=chiapet,D=splicing,E=other),filename=NULL,
                lwd=2,lty=1,col=c(1:5),
                alpha=0.3,cex=3,
                category.names=c("miRNA","tfbs","chiapet","splicing","other"),
                cat.col=c(1:5))
grid.draw(T)

###################################################################################
#Figure 3
###################################################################################

#Following is allele of eQTL overlap with allele of regulation and coding sequence
regulation=read.table("figure/Figure4/data/CAM_alleleRegulatCoding",header=F)
regulation<-regulation[,1]
expression=read.table("figure/Figure4/data/CAM_alleleExpression",header=F)
expression<-expression[,1]
T<-venn.diagram(list(A=regulation,B=expression),filename=NULL,
                lwd=2,col=c("red","blue"),cex=3,lty=1,
                cat.pos=c(-30,30),cat.dist=c(0.05,0.05),cat.cex=2,
                rotation.degree=0
)
grid.draw(T)

#Following is gene of eQTL overlap with gene of regulation and coding sequence
regulation=read.table("figure/Figure4/data/CAM_alleleRegulatCoding",header=F)
regulation<-regulation[,2]
expression=read.table("figure/Figure4/data/CAM_alleleExpression",header=F)
expression<-expression[,2]
T<-venn.diagram(list(A=regulation,B=expression),filename=NULL,
                lwd=2,col=c("red","blue"),cex=3,lty=1,
                cat.pos=c(-30,30),cat.dist=c(0.05,0.05),cat.cex=2,
                rotation.degree=0
)
grid.draw(T)
