#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script find the overlap regions for mutiple bed files and plot a venn diagram(R package eulerr).(Please insure bedtools are in your current PATH)\n',file=stderr())
  cat('Usage:venn_mutipleRegions.R -i=file1,file2,file3 -n=sample1,sample2,sample3 -o=venn.pdf\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\t\tFILE\tThe input bed files separated by comma.(More than 3 and less than 5 files,The files must be sorted by chrom/start.)\n',file=stderr())
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

files=strsplit(inFiles,",")[[1]]
names=strsplit(names,",")[[1]]
if(length(files)==3){
  system(paste("bedtools multiinter -cluster -i ",files[1]," ",files[2]," ",files[3]," >","tmp.multiinter.bed",sep=""))
  multi<-read.delim(file="tmp.multiinter.bed",header=F)
  multi<-multi[,-c(1:5)]
  library("eulerr")
  multi[multi==1]<-"TRUE"
  multi[multi==0]<-"FALSE"
  multi<-matrix(as.logical(as.matrix(multi)),ncol = 3,byrow=F)
  colnames(multi)<-names
  venn<-euler(multi)
  pdf(out)
  plot(venn,counts=T,fill=c("red","green","blue"))
  dev.off()
}else if(length(files)==4){
  system(paste("bedtools multiinter -cluster -i ",files[1]," ",files[2]," ",files[3]," ",files[4]," >","tmp.multiinter.bed",sep=""))
  multi<-read.delim(file="tmp.multiinter.bed",header=F)
  multi<-multi[,-c(1:5)]
  library("eulerr")
  multi[multi==1]<-"TRUE"
  multi[multi==0]<-"FALSE"
  multi<-matrix(as.logical(as.matrix(multi)),ncol = 4,byrow=F)
  colnames(multi)<-names
  venn<-euler(multi)
  pdf(out)
  plot(venn,counts=T,fill=c("red","green","blue","purple"))
  dev.off()
 
}else if(length(files)==5){
  system(paste("bedtools multiinter -cluster -i ",files[1]," ",files[2]," ",files[3]," ",files[4]," ",files[5]," >","tmp.multiinter.bed",sep=""))
  multi<-read.delim(file="tmp.multiinter.bed",header=F)
  multi<-multi[,-c(1:5)]
  library("eulerr")
  multi[multi==1]<-"TRUE"
  multi[multi==0]<-"FALSE"
  multi<-matrix(as.logical(as.matrix(multi)),ncol = 4,byrow=F)
  colnames(multi)<-names
  venn<-euler(multi)
  pdf(out)
  plot(venn,counts=T,fill=c("red","green","blue","purple","yellow"))
  dev.off()
}

system("rm tmp.multiinter.bed")

