#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script plot annotation pie and GO dotpolt for input regions\n',file=stderr())
  cat('Usage:PeakAnnoGO.R peak.bed annotationOut.tsv annoPie.pdf GO.dotplot.pdf\n',file=stderr())
  cat('\t-i\t\tFILE\tInput bed regions\n',file=stderr())
  cat('\t-s\t\tFILE\tThe genome assembly for the input species[mm10,hg19,ce11]\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput prefix for all the output files\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
if(length(args)==0 || args[1]=="-h"){
  usage()
}
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      peakFile=arg.split[2]
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      species=arg.split[2]
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      prefix=arg.split[2]
    }
    if(grepl('^-h',arg)){usage()}
  }
}else{
  usage()
}

if(species == "mm10"){
  library("TxDb.Mmusculus.UCSC.mm10.knownGene",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("org.Mm.eg.db",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
  orgdb<-org.Mm.eg.db
  library("ChIPseeker",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("clusterProfiler",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("ggplot2",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
}else if(species == "hg19"){
  library("TxDb.Hsapiens.UCSC.hg19.knownGene",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("org.Hs.eg.db",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  orgdb<-org.Hs.eg.db
  library("ChIPseeker",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("clusterProfiler",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  library("ggplot2",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
}else if(species == "ce11"){
  library("TxDb.Celegans.UCSC.ce11.refGene",lib.loc="/home/liym/R/lib64/R/library")
  library("org.Ce.eg.db",lib.loc="/home/liym/R/x86_64-pc-linux-gnu-library/3.2")
  txdb<-TxDb.Celegans.UCSC.ce11.refGene
  orgdb<-org.Ce.eg.db
  library("ChIPseeker",lib.loc="/home/liym/R/lib64/R/library")
  library("clusterProfiler",lib.loc="/home/liym/R/lib64/R/library")
  library("ggplot2",lib.loc="/home/liym/R/lib64/R/library")
}else{
  cat('Cannot process the input species',file=stderr())
  q(save="no")
}
peakAnno <- annotatePeak(peakFile,TxDb=txdb)
pdf(file=paste(prefix,".annoPie.pdf",sep=""))
plotAnnoPie(peakAnno)
dev.off()
write.table(as.data.frame(peakAnno),file=paste(prefix,".annotation.tsv",sep=""),quote=F,sep="\t",row.names = F)
gene<-as.data.frame(peakAnno)$geneId
ego <- enrichGO(gene          = gene,
		OrgDb         = orgdb,
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
pdf(file=paste(prefix,".GO.dotplot.pdf",sep=""))
ggplot(subset20,aes(x=p.adjust,y=Description,color=GeneRatio,size=Count))+geom_point(stat='identity')+scale_colour_gradient(high = "#132B43", low = "#56B1F7")
dev.off()
