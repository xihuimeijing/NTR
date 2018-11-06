#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script plot annotation pie and GO dotpolt for input regions\n',file=stderr())
  cat('Usage:PeakAnnoGO.R peak.bed annotationOut.tsv annoPie.pdf GO.dotplot.pdf\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
if(length(args)==0 || args[1]=="-h"){
  usage()
}

library("ChIPseeker")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("clusterProfiler")
library("org.Mm.eg.db")
library("ggplot2")
txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnno <- annotatePeak(args[1],TxDb=txdb)
pdf(file=args[3])
plotAnnoPie(peakAnno)
dev.off()
write.table(as.data.frame(peakAnno),file=args[2],quote=F,sep="\t",row.names = F)
gene<-as.data.frame(peakAnno)$geneId
ego <- enrichGO(gene          = gene,
		OrgDb         = org.Mm.eg.db,
		ont           = "BP",
		pAdjustMethod = "BH",
		pvalueCutoff  = 0.01,
		qvalueCutoff  = 0.05,
		readable      = TRUE)
#Custom dotplot
EGO<-as.data.frame(ego)
EGO$p.adjust<-(-log10(EGO$p.adjust))
ratio<-function(x){as.numeric(strsplit(x,'/')[[1]])[1]/as.numeric(strsplit(x,'/')[[1]])[2]*100}
EGO$GeneRatio<-mapply(FUN=ratio,EGO$GeneRatio,USE.NAMES=F)
subset20<-EGO[c(1:20),]
subset20$Description<-factor(subset20$Description,levels=subset20$Description)
pdf(file=args[4])
ggplot(subset20,aes(x=p.adjust,y=Description,color=GeneRatio,size=Count))+geom_point(stat='identity')+scale_colour_gradient(high = "#132B43", low = "#56B1F7")
dev.off()
