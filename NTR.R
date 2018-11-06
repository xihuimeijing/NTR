#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
usage=function(){
  cat('Usage:Usage:NTR.R -i=H2BGFP.mutiTime.1kb.txt -p=plus -s=1 -e=6 -o=outNTR.txt\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe tab-delimited input file.\n',file=stderr())
  cat('\t-p\tFLOAT\tThe number to be added when the value is zero.\n',file=stderr())
  cat('\t-s\tINT\tThe start column number.[default:1]\n',file=stderr())
  cat('\t-e\tINT\tThe end column number.[default:the last column]\n',file=stderr())
  cat('\t-r\tLOGIC\tThe logical value indicates if need to calculate R-square.[default:false]\n',file=stderr())
  cat('\t-o\tFILE\tThe output file name(Output fields:NTR,pvalue,R-square).\n',file=stderr())
  q(save="no")  
}
start=1
rsquare=FALSE
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }
    if(grepl('^-p=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      plus=as.numeric(arg.split[2])
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      start=as.numeric(arg.split[2])
    }
    if(grepl('^-e=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      end=as.numeric(arg.split[2])
    }
    if(grepl('^-r=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      rsquare=as.logical(arg.split[2])
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
if(! exists("end")){end=ncol(data)}
logData<-log(data[,c(start:end)]+plus)
time=c(0,1,2,4,6,8)
if(rsquare){
  result=matrix(nrow=nrow(data),ncol=3)
  for(i in 1:nrow(data)){
    glm<-glm(as.vector(as.matrix(logData[i,]))~time)
    glmS<-summary(glm)
    result[i,1]=(-1)*glmS$coefficients[2,1]
    result[i,2]=glmS$coefficients[2,4]
    library("vegan")
    result[i,3]= RsquareAdj(glm)$r.squared
  }
}else{
  result=matrix(nrow=nrow(data),ncol=2)
  for(i in 1:nrow(data)){
  	glm<-summary(glm(as.vector(as.matrix(logData[i,]))~time))
  	result[i,1]=(-1)*glm$coefficients[2,1]
  	result[i,2]=glm$coefficients[2,4]
  }
}
result<-round(result,digits=4)
write.table(cbind(data,result),file=out,quote = F, sep ="\t",row.names=F,col.names=F)
