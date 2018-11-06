#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script find the overlap regions for mutiple bed files and plot a venn diagram by using bedr package\n',file=stderr())
  cat('Usage:bedr_venn.R -i=file1,file2,file3 -n=sample1,sample2,sample3 -o=venn.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\t\tFILE\tThe input bed files separated by comma.(More than 3 and less than 5 files,must be sorted by chrom/start)\n',file=stderr())
  cat('\t-n\t\tSTRING\tNames for each files separated by comma.\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput file name[venn.pdf]\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

out="venn.pdf"

if(length(args)>1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -i')
      }else{
        inFiles=arg.split[2]
      }
    }
    if(grepl('^-n=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -n')
      }else{
        names=arg.split[2]
      }
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -o')
      }else{
        out=arg.split[2]
      }
    }
  }
}else if(length(args)==0 || args[1]=="-h"){
  usage()
}
library("bedr")
files=strsplit(inFiles,",")[[1]]
names=strsplit(names,",")[[1]]
peak1<-read.delim(file=files[1],header=F)
peak2<-read.delim(file=files[2],header=F)
peak3<-read.delim(file=files[3],header=F)
filter<-function(peak){
  peak<-peak[,1:3]
  peak$V1<-as.character(peak$V1)
  peak<-peak[!grepl("chrUn|random",peak$V1),]
  return(peak)
}
peak1<-filter(peak1)
peak2<-filter(peak2)
peak3<-filter(peak3)
if(length(files)==3){
  T<-bedr.plot.region(list(names[1]=peak1,names[2]=peak2,names[3]=peak3),filename=NULL,feature="cluster")
  pdf(out)
  grid.draw(T)
  dev.off()
}else if(length(files)==4){
  peak4<-read.delim(file=files[4],header=F)
  peak4<-filter(peak4)
  T<-bedr.plot.region(list(names[1]=peak1,names[2]=peak2,names[3]=peak3,names[4]=peak4),filename=NULL,feature="cluster")
  pdf(out)
  grid.draw(T)
  dev.off()
}else if(length(files)==5){
  peak4<-read.delim(file=files[4],header=F)
  peak4<-filter(peak4)
  peak5<-read.delim(file=files[5],header=F)
  peak5<-filter(peak5)
  T<-bedr.plot.region(list(names[1]=peak1,names[2]=peak2,names[3]=peak3,names[4]=peak4,names[5]=peak5),filename=NULL,feature="cluster")
  pdf(out)
  grid.draw(T)
  dev.off()
}

