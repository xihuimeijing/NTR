#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:Run ChIPseeker package to annotate peaks\n',file=stderr())
  cat('Usage:peakAnnotation.R -i=input.peaks -s=mm10 -h=heatmap.pdf -p=pie.pdf -o=result.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tPeak file in bed/narrowPeak/broadPeak format.\n',file=stderr())
  cat('\t-s\tSTRING\tThe annotation database(hg19/mm10])\n',file=stderr())
  cat('\t-h\tFILE\tHeatmap of peaks to TSS region.[heatmapTSS.pdf]\n',file=stderr())
  cat('\t-p\tFILE\tThe annotation pie file name.[annoPie.pdf]\n',file=stderr())
  cat('\t-o\tFILE\tThe annotation result file if you want to store the result.\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
heatmapFile="heatmapTSS.pdf"
pieFile="annoPie.pdf"
if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the input file -i')
      }else{
        inFile=arg.split[2]
      }
    }else if(grepl('^-s=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the annotation database -s')
      }else{
        database=arg.split[2]
      }
    }else if(grepl('^-h=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      heatmapFile=arg.split[2]
    }else if(grepl('^-p=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      pieFile=arg.split[2]
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      outFile=arg.split[2]
    }
  }
}

library(ChIPseeker)
if(database == "mm10"){
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
}else if(database == "hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
}else{
  stop('Invalid annotation database -s')
}
peakAnno <- annotatePeak(inFile,TxDb=txdb)
pdf(file=pieFile)
plotAnnoPie(peakAnno)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(inFile, windows=promoter)
pdf(file=heatmapFile,width=7,height=14)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()
if(exists("outFile")){
  write.table(as.data.frame(peakAnno),file=outFile,quote=F,sep="\t",row.names = F)
}