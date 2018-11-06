#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Description:This script will sort the input file by means of mutiple columns.\n',file=stderr())
  cat('Usage:sortByMultiColumns.R -i=inputFile.tsv -s=1 -e=6 -o=sortedFile.tsv\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe tab delimted input data to be sorted without header(Can be read from stdin.).\n',file=stderr())
  cat('\t-s\tINT\tThe start column number used to be sorted.\n',file=stderr())
  cat('\t-e\tINT\tThe end column number used to be sorted.\n',file=stderr())
  cat('\t-d\tLOGICAL\tShould the sort order be increasing or decreasing? [default:TRUE for decreasing]\n',file=stderr())
  cat('\t-o\tFILE\tThe sorted output file.\n',file=stderr())
  cat('\t-h\tPrint this help information.\n',file=stderr())
  q(save="no")
}

decending=as.logical("TRUE")
if(length(args)==0 || args[1]=="-h"){
  usage()
}else{
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }else if(grepl('^-s=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the start column number -s')
      }else{
        start=as.numeric(arg.split[2])
      }
    }else if(grepl('^-e=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the start column number -e')
      }else{
        end=as.numeric(arg.split[2])
      }
    }else if(grepl('^-d=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      decending=as.logical(arg.split[2])
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
data<-data[order(rowMeans(data[,c(start:end)]),decreasing = decending),]
write.table(data,file=out,quote=F,sep="\t",row.names = F,col.names = F)
