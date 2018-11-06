#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Usage:rescale_linear.R -i=<file> -s=<INT> -e=<INT> -n=<newrange> -o=<filename>\n',file=stderr())
  cat('Description: This script performs a simple linear conversion of x into the range specified.\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe tab-delimited input file without header.\n',file=stderr())
  cat('\t-s\tINT\tThe start column number to be rescaled\n',file=stderr())
  cat('\t-e\tINT\tThe end column number to be rescaled\n',file=stderr())
  cat('\t-n\tSTRING\tThe newrange separated by comma[default: 0,1]\n',file=stderr())
  cat('\t-o\tFILE\tThe output file name.\n',file=stderr())
  q(save="no")  
}
newRange="0,1"
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      startCol=as.numeric(arg.split[2])
    }
    if(grepl('^-e=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      endCol=as.numeric(arg.split[2])
    }
    if(grepl('^-n=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      newRange=arg.split[2]
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      out=arg.split[2]
    }
    if(grepl('^-h', arg)){
      usage()
    }
  }
}else{
  usage()
}

if(exists("inFile")){
  data<-read.delim(file=inFile,header=F,comment.char='#')
}else{
  data<-read.delim(file('stdin'), header = F)
}
data<-na.omit(data)
new<-as.matrix(data[,c(startCol:endCol)])
minValue=quantile(new,probs=seq(0,1,length.out=101))[2]
maxValue=quantile(new,probs=seq(0,1,length.out=101))[99]
new[new>maxValue]=maxValue
new[new<minValue]=minValue
library(plotrix)
newStart<-as.numeric(strsplit(newRange,",")[[1]][1])
newEnd<-as.numeric(strsplit(newRange,",")[[1]][2])
new<-rescale(new,newrange=c(newStart,newEnd))
new<-cbind(data[,c(1:(startCol-1))],new)
write.table(new,file=out,sep="\t",quote=F,row.names=F,col.names=F)

