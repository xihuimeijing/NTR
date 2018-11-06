#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script will normalize each column by the mean of this column\n',file=stderr())
  cat('Usage:normalizedByColMeans.R -i=inputFile.tsv -o=normalizedFile.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe tab delimted input data without header(Can be read from stdin.).\n',file=stderr())
  cat('\t-s\tINT\tThe start column number.[default:1]\n',file=stderr())
  cat('\t-e\tINT\tThe end column number.[default:the last column]\n',file=stderr())
  cat('\t-n\tSTRING\tThe normalization method(subtract,ratio).[default:ratio]\n',file=stderr())
  cat('\t-o\tFILE\tThe normalized output file.\n',file=stderr())
  cat('\t-h\tPrint this help information.\n',file=stderr())
  q(save="no")
}
start=1
norm="ratio"
if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }else if(grepl('^-s=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        start=as.numeric(arg.split[2])
      }else if(grepl('^-e=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        ed=as.numeric(arg.split[2])
      }else if(grepl('^-n=', arg)){
        arg.split = strsplit(arg, '=', fixed = T)[[1]]
        norm=arg.split[2]
    }else if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the output file name -o')
      }else{
        out=arg.split[2]
      }
    }
  }
}

if(exists("inFile")){
    data<-read.delim(file=inFile,header=F,comment.char='#')
}else{
    data<-read.delim(file('stdin'), header = F)
}

if(! exists("ed")){ed=ncol(data)}
scaler<-colMeans(data[,c(start:ed)])
if(norm=="ratio"){
  result<-data[,c(start:ed)]/t(replicate(nrow(data),scaler))
}
if(norm=="subtract"){
  result<-data[,c(start:ed)]-t(replicate(nrow(data),scaler))
}
write.table(cbind(data,result),file=out,row.names = F,col.names = F,quote=F,sep="\t")