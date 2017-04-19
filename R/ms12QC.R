
##' @title Calculate the MS1 and MS2 level QC metrics
##' @description Calculate the MS1 level QC metrics
##' @param spectraList An experiment design input file
##' @param outdir Output directory
##' @param cpu The number of cpu used
##' @return A data frame 
##' @export
##' @author Bo Wen \email{wenbo@@genomics.cn}
calcMSQCMetrics=function(spectraList=NULL,cpu=2,outdir="./"){
  
  #library("Rcpp")
  exp <- read.delim(spectraList)
  if(cpu==0){
    cpu <- detectCores()
    
  }
  if(cpu>=4){
    cpu=4
  }
  cl <- makeCluster(getOption("cl.cores", cpu))
  clusterEvalQ(cl,library("MSnbase"))
  clusterExport(cl, c("outdir"),envir=environment())
  res<-parSapply(cl,as.character(exp[,.INPUT.COLS["FILE"]]),function(x){
    mz <- readMSData(x,msLevel=1)
    outfile <- paste(outdir,"/",basename(x),"-ms1qc.txt",collapse="",sep="")
    write.table(x=header(mz),file=outfile,quote=FALSE,sep="\t",
                col.names=TRUE, row.names=FALSE)
    
    mz2 <- readMSData(x,msLevel=2)
    outfile2 <- paste(outdir,"/",basename(x),"-ms2qc.txt",collapse="",sep="")
    write.table(x=header(mz2),file=outfile2,quote=FALSE,sep="\t",
                col.names=TRUE, row.names=FALSE)
    out <- c(outfile,outfile2)
    ## return a matrix. Can not use data.frame here.
    out
    
  })
  stopCluster(cl)
  res<-as.data.frame(t(res))
  names(res) <- c("MS1QC","MS2QC")
  res[,"file"]=rownames(res)
  res <- merge(res,exp,by=.INPUT.COLS["FILE"])
  
  res
  
}

plotMS1TIC=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   y = read.delim(MS1QC)$tic,
                   x= read.delim(MS1QC)$retention.time/60)
  
  ggplot.RT(p,fig=fig,xlab="Retention time",ylab="TIC")
}

plotMS1PeaksCount=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   y = read.delim(MS1QC)$peaks.count,
                   x= read.delim(MS1QC)$retention.time/60)
  
  ggplot.RT(p,fig=fig,xlab="Retention time",ylab="Peaks count")
}

plotMS1IonCount=function(x,fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   y = read.delim(MS1QC)$ionCount,
                   x= read.delim(MS1QC)$retention.time/60)
  
  ggplot.RT(p,fig=fig,xlab="Retention time",ylab="Ion count")
}

plotMS2PeakFreq=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   x= read.delim(MS2QC)$retention.time/60)
  png(fig,width=1200,height=1000,res=200)
  gg.obj <- ggplot(data=p,aes(x=x,colour=as.factor(techRep),
                                 linetype=as.factor(bioRep))) +
    geom_density(size=0.2,alpha=0.8)+
    #scale_x_continuous(labels = comma)+
    #scale_y_continuous(labels=comma)+
    theme(axis.text.x  = element_text(size=6,angle=90,vjust=0.5))+
    theme(axis.text.y  = element_text(size=6))+
    xlab(label="Retention time")+
    #ylab(label="rt")+
    facet_wrap( ~ fraction+sample, ncol = 6)+
    labs(colour="Technical replicate",linetype="Biological replicate")+
    coord_flip()+
    theme(legend.text=element_text(size=6.5),legend.title=element_text(size=7))
  
  print(gg.obj)
  dev.off()
  

}

plotMS1boxplot=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   x= read.delim(MS1QC)$retention.time/60)
  png(fig,width=1200,height=1000,res=200)
  gg.obj <- ggplot(data=p,aes(y=x,x=as.factor(fraction))) +
    geom_boxplot(size=0.2,alpha=0.8)+
    #scale_x_continuous(labels = comma)+
    #scale_y_continuous(labels=comma)+
    theme(axis.text.x  = element_text(size=6,angle=90,vjust=0.5))+
    theme(axis.text.y  = element_text(size=6))+
    ylab(label="Retention time")+
    xlab(label="Fraction")+
    facet_wrap( ~ sample+bioRep+techRep, ncol = 3)+
    #labs(colour="Technical replicate",linetype="Biological replicate")+
    #coord_flip()+
    theme(legend.text=element_text(size=6.5),legend.title=element_text(size=7))
  
  print(gg.obj)
  dev.off()
}

plotMS1Count=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  pdat <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   ms1count= nrow(read.delim(MS1QC)))
  pdat$sample <- as.character(pdat$sample)
  pdat$bioRep <- as.character(pdat$bioRep)
  pdat$techRep <- as.character(pdat$techRep)
  #pdat$sample <- as.character(pdat$sample)
  plotClass <- "ms1count"
  x<-reshape2::dcast(pdat,
                     sample+bioRep+fraction~techRep,
                     fill=0,value.var=c(plotClass))
  
  m<-reshape2::melt(x,id.vars=c("sample","bioRep","fraction"),
                    value.name=c(plotClass),
                    variable.name="techRep")
  
  #if(length(unique(m$techRep))>=6){
  if(max(nchar(levels(as.factor(m$techRep))))>=6){
    rotate=90 
  }else{ 
    rotate=0
  }
  
  
  png(fig,width=1000,height=800,res=200)
  
  p<-ggplot(m,aes_string(x="fraction",y=plotClass,
                         linetype="techRep",
                         colour="sample",shape="bioRep"))+
    geom_point()+
    geom_line()+
    ylab("MS1 Count")+
    xlab("Fraction")+
    expand_limits(y=0)
  print(p)
  
  dev.off()
  #pngFile <- paste(basename(res$input_parameter$report.dir),basename(pngFile),
  #                 sep="/")
  #return(fig)
  
  
}


plotMS1CountErrorBar=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  pdat <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                      summarise,
                      ms1count= nrow(read.delim(MS1QC)))
  pdat$sample <- as.character(pdat$sample)
  pdat$bioRep <- as.character(pdat$bioRep)
  pdat$techRep <- as.character(pdat$techRep)
  #pdat$sample <- as.character(pdat$sample)
  plotClass <- "ms1count"
  
  z<-ddply(pdat,.(sample,fraction),
           function(x){
             data.frame(val=mean(x[,plotClass]),
                        se=sd(x[,plotClass])/sqrt(length(x[,plotClass])))})
  
  #pngFile=paste(outdir,"/","fig_ms1_sample_error_bar_",plotClass,".png",
  #              sep="",collapse="")
  png(fig,width=1000,height=800,res=200)
  p<-ggplot(z,aes(x=fraction, y = val,fill=sample,colour=sample))+
    expand_limits(y = 0)+
    geom_line()+
    geom_point()+
    ylab("MS1 Count")+
    xlab("Fraction")+
    geom_errorbar(aes(ymin=val-se, ymax=val+se),width=0.3)+
    #scale_x_continuous(breaks = seq(1,max(res_fraction_level$fraction),1))+
    #theme(axis.text.x  = element_text(angle=0))+
    #axis.title=element_text(face="bold",size=15),
    #plot.title=element_text(face="bold",size=20))+
    scale_fill_hue(c=90,l=50)
  #if(res$input_parameter$maxFraction > 6){
   # p <- p + theme(axis.text.x  = element_text(angle=90,vjust=0.5 ))
  #}
  print(p)
  dev.off()
  #pngFile <- paste(basename(res$input_parameter$report.dir),basename(pngFile),
  #                 sep="/")
  return(fig)
  
}


plotMS2boxplot=function(x, fig="test.png"){
  x$MS1QC <- as.character(x$MS1QC)
  x$MS2QC <- as.character(x$MS2QC)
  p <- plyr::ddply(x,.(sample,bioRep,techRep,fraction),
                   summarise,
                   x= read.delim(MS2QC)$retention.time/60)
  png(fig,width=1200,height=1000,res=200)
  gg.obj <- ggplot(data=p,aes(y=x,x=as.factor(fraction))) +
    geom_boxplot(size=0.2,alpha=0.8)+
    #scale_x_continuous(labels = comma)+
    #scale_y_continuous(labels=comma)+
    theme(axis.text.x  = element_text(size=6,angle=90,vjust=0.5))+
    theme(axis.text.y  = element_text(size=6))+
    ylab(label="Retention time")+
    xlab(label="Fraction")+
    facet_wrap( ~ sample+bioRep+techRep, ncol = 3)+
    #labs(colour="Technical replicate",linetype="Biological replicate")+
    #coord_flip()+
    theme(legend.text=element_text(size=6.5),legend.title=element_text(size=7))
  
  print(gg.obj)
  dev.off()
}



ggplot.RT=function(data=NULL,fig=NULL,xlab=NULL,ylab=NULL){
  png(fig,width=1200,height=1000,res=200)
  gg.obj <- ggplot(data=data,aes(x=x,y=y,colour=as.factor(techRep),
                              linetype=as.factor(bioRep))) +
    geom_line(size=0.2,alpha=0.8)+
    #scale_x_continuous(labels = comma)+
    #scale_y_continuous(labels=comma)+
    theme(axis.text.x  = element_text(size=6,angle=90,vjust=0.5))+
    theme(axis.text.y  = element_text(size=6))+
    xlab(label=xlab)+
    ylab(label=ylab)+
    #facet_wrap( ~ fraction, ncol = 6)+
    facet_wrap( ~ fraction+sample, ncol = 6)+
    labs(colour="Technical replicate",linetype="Biological replicate")+
    coord_flip()+
    theme(legend.text=element_text(size=6.5),legend.title=element_text(size=7))
    
  
  
  print(gg.obj)
  dev.off()
  
}

plotMS12=function(res=NULL, outdir="./"){
  outfig <- list()
  outfig$ms1tic <- paste(outdir,"/ms1tic.png",sep="")
  outfig$ms1peakscount <- paste(outdir,"/ms1peakscount.png",sep="")
  outfig$ms1ioncount <- paste(outdir,"/ms1ioncount.png",sep="")
  outfig$ms2peaksdensity <- paste(outdir,"/ms2peaksdensity.png",sep="")
  outfig$ms2boxplot <- paste(outdir,"/ms2boxplot.png",sep="")
  outfig$ms1countdot <- paste(outdir,"/ms1countdot.png",sep="")
  outfig$ms1counterrorbar <- paste(outdir,"/ms1counterrorbar.png",sep="")
  outfig$ms1boxplot <- paste(outdir,"/ms1boxplot.png",sep="")
  
  plotMS1TIC(res,fig=outfig$ms1tic)
  plotMS1PeaksCount(res,fig=outfig$ms1peakscount)
  plotMS1IonCount(res,fig=outfig$ms1ioncount)
  plotMS2PeakFreq(res,fig=outfig$ms2peaksdensity)
  plotMS2boxplot(res,fig=outfig$ms2boxplot)
  plotMS1boxplot(res,fig=outfig$ms1boxplot)
  plotMS1Count(res,fig=outfig$ms1countdot)
  plotMS1CountErrorBar(res,fig=outfig$ms1counterrorbar)
  
  return(outfig)
}
