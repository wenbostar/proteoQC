
##' @title Barplot in different level for each fraction
##' @description Barplot in different level for each fraction
##' @param res An object of msQCres
##' @param level 1: total spectrum, 2: identified spectrum, 
##' 3: identified peptide, 4: identified protein.
##' @return The name of the figure
plotFractionIDResult=function(res,level=NA){
  
  
  res_fraction_level <- res$res_fraction_level
  outdir <- res$input_parameter$report.dir
  plotClass<- switch(level,
         "1"=.RES.FRACTION.COLS["SPECTRUM_TOTAL"],
         "2"=.RES.FRACTION.COLS["SPECTRUM"],
         "3"=.RES.FRACTION.COLS["PEPTIDE"],
         "4"=.RES.FRACTION.COLS["PROTEIN"],
         stop("The value of level should be 1,2,3"))
  
  
  res_fraction_level[,.INPUT.COLS['SAMPLE']]<-
    as.character(res_fraction_level[,.INPUT.COLS['SAMPLE']])
  res_fraction_level[,.INPUT.COLS['BIOREP']]<-
    as.character(res_fraction_level[,.INPUT.COLS['BIOREP']])
  res_fraction_level[,.INPUT.COLS['TECHREP']]<-
    as.character(res_fraction_level[,.INPUT.COLS['TECHREP']])
  colnames(res_fraction_level[,c(.INPUT.COLS['SAMPLE'],
                                 .INPUT.COLS['BIOREP'],
                                 .INPUT.COLS['TECHREP'],
                                 .INPUT.COLS['FRACTION'],
                                 .RES.FRACTION.COLS['SPECTRUM_TOTAL'],
                                 .RES.FRACTION.COLS['SPECTRUM'],
                                 .RES.FRACTION.COLS['PEPTIDE'],
                                 .RES.FRACTION.COLS['PROTEIN']
                                 )])<-c("sample",
                                        "bioRep",
                                        "techRep",
                                        "fraction",
                                        "spectrum_total",
                                        "spectrum",
                                        "peptide",
                                        "protein")
  

  x<-reshape2::dcast(res_fraction_level,
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
  
  pngFile=paste(outdir,"/","fig_fraction_bar_",plotClass,".png",
                sep="",collapse="")
  png(pngFile,width=1000,height=800,res=200)
  p<-ggplot(m,aes_string(x="fraction",y=plotClass,
                         linetype="techRep",
                         colour="sample",shape="bioRep"))+
    geom_point()+
    geom_line()+
    expand_limits(y=0)
  print(p)
  
  dev.off()
  pngFile <- paste(basename(res$input_parameter$report.dir),basename(pngFile),
                   sep="/")
  return(pngFile)

}

##' @title Error barplot in different level for each fraction
##' @description Error Barplot in different level for each fraction
##' @param res An object of parser result
##' @param level 1: total spectrum, 2: identified spectrum, 
##' 3: identified peptide, 4: identified protein.
##' @return The name of the figure
plotSampleIDResultErrorBar=function(res,level=NA){
  
  res_fraction_level <- res$res_fraction_level
  outdir <- res$input_parameter$report.dir
  plotClass<- switch(level,
                     "1"=.RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                     "2"=.RES.FRACTION.COLS["SPECTRUM"],
                     "3"=.RES.FRACTION.COLS["PEPTIDE"],
                     "4"=.RES.FRACTION.COLS["PROTEIN"],
                     stop("The value of level should be 1,2,3,4"))
  
  res_fraction_level[,.INPUT.COLS['SAMPLE']]<-
    as.character(res_fraction_level[,.INPUT.COLS['SAMPLE']])
  res_fraction_level[,.INPUT.COLS['BIOREP']]<-
    as.character(res_fraction_level[,.INPUT.COLS['BIOREP']])
  res_fraction_level[,.INPUT.COLS['TECHREP']]<-
    as.character(res_fraction_level[,.INPUT.COLS['TECHREP']])
  colnames(res_fraction_level[,c(.INPUT.COLS['SAMPLE'],
                                 .INPUT.COLS['BIOREP'],
                                 .INPUT.COLS['TECHREP'],
                                 .INPUT.COLS['FRACTION'],
                                 .RES.FRACTION.COLS['SPECTRUM_TOTAL'],
                                 .RES.FRACTION.COLS['SPECTRUM'],
                                 .RES.FRACTION.COLS['PEPTIDE'],
                                 .RES.FRACTION.COLS['PROTEIN']
  )])<-c("sample",
         "bioRep",
         "techRep",
         "fraction",
         "spectrum_total",
         "spectrum",
         "peptide",
         "protein")
  

  z<-ddply(res_fraction_level,.(sample,fraction),
           function(x){
             data.frame(val=mean(x[,plotClass]),
                        se=sd(x[,plotClass])/sqrt(length(x[,plotClass])))})
  
  pngFile=paste(outdir,"/","fig_sample_error_bar_",plotClass,".png",
                sep="",collapse="")
  png(pngFile,width=1000,height=800,res=200)
  y.lab <- switch(level,
                     "1"="Total spectra",
                     "2"="Identified spectra",
                     "3"="Identified peptides",
                     "4"="Identified proteins",
                     stop("The value of level should be 1,2,3,4"))
  p<-ggplot(z,aes(x=fraction, y = val,fill=sample,colour=sample))+
    expand_limits(y = 0)+
    geom_line()+
    geom_point()+
    ylab(y.lab)+
    xlab("Fraction")+
    geom_errorbar(aes(ymin=val-se, ymax=val+se),width=0.3)+
    #scale_x_continuous(breaks = seq(1,max(res_fraction_level$fraction),1))+
    #theme(axis.text.x  = element_text(angle=0))+
                     #axis.title=element_text(face="bold",size=15),
                     #plot.title=element_text(face="bold",size=20))+
    scale_fill_hue(c=90,l=50)
  if(res$input_parameter$maxFraction > 6){
    p <- p + theme(axis.text.x  = element_text(angle=90,vjust=0.5 ))
  }
  print(p)
  dev.off()
  pngFile <- paste(basename(res$input_parameter$report.dir),basename(pngFile),
                   sep="/")
  return(pngFile)
  
}

##' @title plot MS2 mass error
##' @description plot MS2 mass error
##' @param res An object of msQCres
##' @return The name of the figure
plotMS2Error_obsolete <- function(res) {
  outdir <- res$input_parameter$report.dir
  #load("result.rda")    
  p <- plyr::ddply(res$res_fraction_level,
                   .(sample,bioRep,techRep,fraction),
                   summarise,
                   error = read.delim(peptide_summary, 
                                      stringsAsFactor=FALSE,
                                      colClasses=.PepSummaryColClass)$ms2delta)
  f <- plyr::ddply(res$res_fraction_level,
                   .(sample,bioRep,techRep,fraction),
                   function(x) 
                   quantile(as.numeric(read.delim(x[,"peptide_summary"], 
                                      stringsAsFactor=FALSE,
                                      colClasses=.PepSummaryColClass)$ms2delta),
                            probs = seq(0, 1, 0.05))[c("5%","95%")])

  fig <- paste(outdir,"/","ms2_error.png",sep="",collapse="")
  png(fig,width=1000,height=800,res=200)
  qcHist(p,f,xlab="Fragment ion mass error (Da)")
  dev.off()
  fig <- paste(basename(res$input_parameter$report.dir),
               basename(fig),
               sep = "/")
  return(fig)  
} 


##' @title plot MS2 mass error
##' @description plot MS2 mass error
##' @param res An object of msQCres
##' @return The name of the figure
plotMS2Error <- function(res) {
  outdir <- res$input_parameter$report.dir
  #load("result.rda")    
  p <- plyr::ddply(res$res_fraction_level,
                   .(sample,bioRep,techRep,fraction),
                   function(f){
                   stats::fivenum(as.numeric(unlist(strsplit(
                     x=as.character(read.delim(f$peptide_summary)$ms2delta),
                     split=";"))))
                   })
  p$techRep <- as.character(p$techRep)
  p$sample <- as.character(p$sample)
  p$bioRep <- as.character(p$bioRep)
  
  fig <- paste(outdir,"/","ms2_error.png",sep="",collapse="")
  png(fig,width=1000,height=800,res=200)
  
  gg <- ggplot(data=p,aes(x=fraction,
                          linetype=techRep,
                          colour=sample,shape=bioRep))+
    geom_point(aes(y=V1))+geom_line(aes(y=V1))+
    geom_point(aes(y=V2))+geom_line(aes(y=V2))+
    geom_point(aes(y=V3))+geom_line(aes(y=V3))+
    geom_point(aes(y=V4))+geom_line(aes(y=V4))+
    geom_point(aes(y=V5))+geom_line(aes(y=V5))+
    ylab("Fragment ion mass error(Da)")
  print(gg)
  dev.off()
  fig <- paste(basename(res$input_parameter$report.dir),
               basename(fig),
               sep = "/")
  return(fig)  
} 


##' @title plot MS1 mass error
##' @description plot MS1 mass error
##' @param res An object of msQCres
##' @param plot.class ppm or da
##' @return The name of the figure
plotMS1Error=function(res,plot.class="ppm") {  
    outdir <- res$input_parameter$report.dir    
    err_up=0.05
    err_down=-0.05
    step=30 
    br=seq(err_down,err_up,(err_up-err_down)/step)
    curenv <<- environment()  
    if (plot.class == "ppm") {
        p<-ddply(res$res_fraction_level,
                 .(sample,bioRep,techRep,fraction),
                 summarise,
                 error = {
                     peps <- read.delim(peptide_summary, 
                                        stringsAsFactors = FALSE,
                                        colClasses=.PepSummaryColClass);
                     delta.da  <- as.numeric(peps$delta_da)
                     delta.ppm <- as.numeric(peps$delta_ppm)
                     restmp <- eval( get("res" , envir=curenv ))
                     if (restmp$input_parameter$tolu == "ppm") {    
                         delta.da  <- delta.da[abs(delta.ppm) <= 
                                                restmp$input_parameter$tol]
                         delta.ppm <- delta.ppm[abs(delta.ppm) <= 
                                                restmp$input_parameter$tol] 
                     } else {
                         ## dalton
                         delta.da  <- delta.da[abs(delta.da) <= 
                                                 restmp$input_parameter$tol]
                         delta.ppm <- delta.ppm[abs(delta.da) <= 
                                                  restmp$input_parameter$tol]                     
                     }
                     delta.ppm
                 })
        f <- ddply(res$res_fraction_level,
                   .(sample,bioRep,techRep,fraction), function(x) {
                       peps <-read.delim(x[1,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                                          colClasses=.PepSummaryColClass,
                                          stringsAsFactors = FALSE);
                       ##peps <- read.delim(peptide_summary,stringsAsFactors=F);
                       delta.da  <- as.numeric(peps$delta_da)
                       delta.ppm <- as.numeric(peps$delta_ppm)
                       restmp=eval( get("res" , envir=curenv ))
                       if (restmp$input_parameter$tolu == "ppm") {        
                           delta.da  <- delta.da[abs(delta.ppm) <= 
                                                   restmp$input_parameter$tol]
                           delta.ppm <- delta.ppm[abs(delta.ppm) <= 
                                                    restmp$input_parameter$tol]        
                       } else {
                           ## dalton
                           delta.da  <- delta.da[abs(delta.da) <= 
                                                   restmp$input_parameter$tol]
                           delta.ppm <- delta.ppm[abs(delta.da) <= 
                                                    restmp$input_parameter$tol]        
                       }
                       quantile(delta.ppm, 
                                probs = seq(0, 1, 0.05))[c("5%","95%")]
                   })
        fig <- paste(outdir,"/","ms1_error_ppm.png",sep="",collapse="")
        png(fig,width=1000,height=800,res=200)
        qcHist(p,f,xlab="Precusor ion mass error (ppm)")
        dev.off()
        fig <- paste(basename(res$input_parameter$report.dir),basename(fig),
                     sep="/")
        return(fig)
    } else {
        p<-ddply(res$res_fraction_level,.(sample,bioRep,techRep,fraction),
                 summarise,
                 error = {
                     peps <- read.delim(peptide_summary, 
                                        stringsAsFactors = FALSE,
                                        colClasses=.PepSummaryColClass);
                     delta.da  <- as.numeric(peps$delta_da)
                     delta.ppm <- as.numeric(peps$delta_ppm)
                     restmp=eval( get("res" , envir=curenv ))
                     if (restmp$input_parameter$tolu == "ppm"){                     
                         delta.da  <- delta.da[abs(delta.ppm) <= 
                                                 restmp$input_parameter$tol]
                         delta.ppm <- delta.ppm[abs(delta.ppm) <= 
                                                  restmp$input_parameter$tol]                     
                     } else {
                         ## dalton
                         delta.da  <- delta.da[abs(delta.da) <= 
                                                 restmp$input_parameter$tol]
                         delta.ppm <- delta.ppm[abs(delta.da) <= 
                                                  restmp$input_parameter$tol]
                         
                     }
                     delta.da
                 })
        f<-ddply(res$res_fraction_level,.(sample,bioRep,techRep,fraction),
                 function(x){
            peps <- read.delim(x[1,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                               stringsAsFactors=FALSE,
                               colClasses=.PepSummaryColClass);
            ## peps <- read.delim(peptide_summary,stringsAsFactors=F);
            delta.da  <- as.numeric(peps$delta_da)
            delta.ppm <- as.numeric(peps$delta_ppm)
            restmp=eval( get("res" , envir=curenv ))
            if (restmp$input_parameter$tolu == "ppm") {            
                delta.da  <-delta.da[abs(delta.ppm)<=restmp$input_parameter$tol]
                delta.ppm<-delta.ppm[abs(delta.ppm)<=restmp$input_parameter$tol]            
            } else {
                ## dalton
                delta.da  <- delta.da[abs(delta.da)<=restmp$input_parameter$tol]
                delta.ppm <-delta.ppm[abs(delta.da)<=restmp$input_parameter$tol]            
            }
            quantile(delta.da, probs = seq(0, 1, 0.05))[c("5%","95%")]
        })
        ## save(p,f,file="pf.rda")
        fig <- paste(outdir,"/","ms1_error_da.png",sep="",collapse="")
        png(fig,width=1000,height=800,res=200)
        qcHist(p,f,xlab="Precusor ion mass error (Da)")
        dev.off()
        fig <- paste(basename(res$input_parameter$report.dir), basename(fig), 
                     sep="/")
        return(fig)
    }  
}

qcHist2<-function(p,f,xlab=NA){
  save(p,f,xlab,file="ok.rda")
  m<-melt(f,id.vars=c("sample","bioRep","techRep","fraction"),
          value.name="fractile",
          variable.name="label")

  if(length(unique(p$fraction))>12){
    gp=ggplot(p,aes(x=error))+
      geom_histogram()+
      geom_vline(data=m,aes(xintercept = fractile),colour="red", 
                 linetype = "longdash")+
      facet_wrap(sample+bioRep+techRep~fraction,ncol=12)+
      theme_bw()+
      theme(axis.text.x  = element_text(angle=90,vjust=0.5),
            axis.title=element_text(face="bold",size=15),
            plot.title=element_text(face="bold",size=20))
    max_y=max(ggplot_build(gp)$data[[1]]$count)  
    gp <- gp+scale_y_continuous(expand = c(0,0),limits=c(0,max_y*1.2))+
      geom_text(data=m,aes(x=fractile,label=round(fractile,digits=3)),
                y=max_y*1.1,size=4,colour="blue")
    if(!is.na(xlab)){
      gp <- gp+ggplot2::xlab(label=xlab)
    }
    print(gp)
    
  }else{

    gp=ggplot(p,aes(x=error))+
      geom_histogram()+
      geom_vline(data=m,aes(xintercept = fractile),colour="red", 
                 linetype = "longdash")+
      facet_grid(sample+bioRep+techRep~fraction)+
      theme_bw()+
      theme(axis.text.x  = element_text(angle=90,vjust=0.5),
            axis.title=element_text(face="bold",size=15),
            plot.title=element_text(face="bold",size=20))

    max_y=max(ggplot_build(gp)$data[[1]]$count)	
  
    gp <- gp+scale_y_continuous(limits=c(0,max_y*1.3))+
      geom_text(data=m,aes(x=fractile,label=round(fractile,digits=3)),
                y=max_y*1.1,size=4,colour="blue")
    if(!is.na(xlab)){
      gp <- gp+ggplot2::xlab(label=xlab)
    }
    print(gp)
    
  }
}


qcHist<-function(p,f,xlab=NA){
  p$sample <- as.character(p$sample)
  p$bioRep <- as.character(p$bioRep)
  p$techRep <- as.character(p$techRep)
  
  gg.obj <- ggplot(data=p,aes(x=error,linetype=techRep,
                              size=bioRep,colour=sample)) +
    #geom_density(size=0.6)+
    geom_density()+
    geom_vline(aes(xintercept=0),alpha=0.3,size=1.2)+
    theme(axis.text.x  = element_text(size=6,angle=90,vjust=0.5))+
    theme(axis.text.y  = element_text(size=6))+
    facet_wrap( ~ fraction, ncol = 6)+
    labs(linetype="Technical replicate",size="Biological replicate",
         colour="Sample")+
    coord_flip()#+
    #theme(legend.text=element_text(size=0.7),legend.title=element_text(size=7))
  if(!is.na(xlab)){
    gg.obj <- gg.obj + xlab(label=xlab)
  }
  print(gg.obj)
}


##' @title Venn plot in technical replicate level
##' @description Venn plot in technical replicate level
##' @param res An object of msQCres
##' @return The name of the figure
plotTechRepVenn=function(res){
    
  outdir <- res$input_parameter$report.dir
  result <- ddply(res$res_techRep_level,
                  .variables=c(.INPUT.COLS["SAMPLE"],
                               .INPUT.COLS["BIOREP"]),function(x){
    ## when the number of technical replicates are >=2, it will work.
    if(dim(x)[1]>=2){
      ## peptide level
      venn.data <-dlply(x,.variables=c(.INPUT.COLS["TECHREP"]),function(y){
        pep <- read.delim(y[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                          stringsAsFactors=FALSE,
                          colClasses=.PepSummaryColClass)
        pep <- unique(pep[,'peptide'])
        
      })
      
      #vobj <- Venn(venn.data)
      
      #par(bty="n")
      
      ven.fig.pep <- paste(outdir,"/",x[1,.INPUT.COLS["SAMPLE"]],"-",
                           x[1,.INPUT.COLS["BIOREP"]],"-pepvenn.png",
                           sep="",collapse="")
      png(ven.fig.pep,width=800,height=800,res=200)
      par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
      #plot(vobj, doWeights = FALSE)
      #venn.diagram(x=venn.data,filename=ven.fig.pep,
      #cex=1.1,cat.cex=1.4,reverse=TRUE,resolution=200,height=700,
      #width=700,compression="none")
      grob.list <- venn.diagram(x=venn.data,filename=NULL,
                                cex=1.1,cat.cex=1.4,reverse=TRUE,
                                resolution=200,height=700,width=700,
                                compression="none")
      grid.draw(grob.list)
      dev.off()
      ven.fig.pep <- paste(basename(outdir),basename(ven.fig.pep),sep="/")
      
      ## 2. protein level
      venn2.data <-dlply(x,.variables=c(.INPUT.COLS["TECHREP"]),function(y){
        pro <- read.delim(y[,.RES.FRACTION.COLS["PRO_SUMMARY"]],
                          colClasses=.ProSummaryColClass)
        pro <- unique(pro[,'Accession'])
        
      })
      
      #vobj2 <- Venn(venn2.data)
      #par(bty="n")
      ven.fig.pro <- paste(outdir,"/",x[1,.INPUT.COLS["SAMPLE"]],"-",
                           x[1,.INPUT.COLS["BIOREP"]],"-provenn.png",
                           sep="",collapse="")
      png(ven.fig.pro,width=800,height=800,res=200)
      #plot(vobj2, doWeights = FALSE)
      par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
      grob.list <- venn.diagram(x=venn2.data,filename=NULL,
                                cex=1.1,cat.cex=1.4,reverse=TRUE,
                                resolution=200,height=700,width=700,
                                compression="none")
      grid.draw(grob.list)
      dev.off()
      ven.fig.pro <- paste(basename(outdir),basename(ven.fig.pro),sep="/")
      data.frame(pepfig=ven.fig.pep,profig=ven.fig.pro)
    }else{
      data.frame(pepfig=NA,profig=NA)
    }
    
  })

  return(result)
  
}



##' @title Venn plot in biological replicate level
##' @description Venn plot in biological replicate level
##' @param res An object of msQCres
##' @return The name of the figure
plotBioRepVenn=function(res){
  
  outdir <- res$input_parameter$report.dir
  result <- ddply(res$res_bioRep_level,
                  .variables=c(.INPUT.COLS["SAMPLE"]),function(x){
    # when the number of biological replicates are >=2, it will work.
    if(dim(x)[1]>=2){
      ## 1. peptide level
      venn.data <-dlply(x,.variables=c(.INPUT.COLS["BIOREP"]),function(y){
        pep <- read.delim(y[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                          colClasses=.PepSummaryColClass)
        pep <- unique(pep[,'peptide'])
        
      })
      

      #par(bty="n")
      
      ven.fig.pep <- paste(outdir,"/",x[1,.INPUT.COLS["SAMPLE"]],"-pepvenn.png",
                           sep="",collapse="")
      png(ven.fig.pep,width=800,height=800,res=200)
      #plot(vobj, doWeights = FALSE)
      par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
      grob.list <- venn.diagram(x=venn.data,filename=NULL,
                                cex=1.1,cat.cex=1.4,reverse=TRUE,
                                resolution=200,height=700,width=700,
                                compression="none")
      grid.draw(grob.list)
      dev.off()
      ven.fig.pep <- paste(basename(outdir),basename(ven.fig.pep),sep="/")
      
      ## 2. protein level
      venn2.data <-dlply(x,.variables=c(.INPUT.COLS["BIOREP"]),function(y){
        pro <- read.delim(y[,.RES.FRACTION.COLS["PRO_SUMMARY"]],
                          colClasses=.ProSummaryColClass)
        pro <- unique(pro[,'Accession'])
        
      })
      
      #vobj2 <- Venn(venn2.data)
      #par(bty="n")
      ven.fig.pro <- paste(outdir,"/",x[1,.INPUT.COLS["SAMPLE"]],"-provenn.png",
                           sep="",collapse="")
      png(ven.fig.pro,width=500,height=500,res=100)
      par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
      grob.list <- venn.diagram(x=venn2.data,filename=NULL,cex=1.1,
                                cat.cex=1.4,reverse=TRUE,resolution=200,
                                height=700,width=700,compression="none")
      grid.draw(grob.list)
      #plot(vobj2, doWeights = FALSE)
      dev.off()
      ven.fig.pro <- paste(basename(outdir),basename(ven.fig.pro),sep="/")
      data.frame(pepfig=ven.fig.pep,profig=ven.fig.pro)
    }else{
      data.frame(pepfig=NA,profig=NA)
    }
    
  })

  return(result)

}

##' @title Venn plot in sample level
##' @description Venn plot in sample level
##' @param res An object of msQCres
##' @return The name of the figure
plotSampleVenn = function(res) {
  outdir <- res$input_parameter$report.dir
  result <- data.frame()
  if (res$input_parameter$maxSample>=2) {
    # when the number of samples are >=2, it will work.
    ## 1. peptide level
    venn.data <- dlply(res$res_sample_level,
                       .variables=c(.INPUT.COLS["SAMPLE"]),function(y){
      pep <- read.delim(y[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                        colClasses=.PepSummaryColClass)
      pep <- unique(pep[,'peptide'])      
    })
    
    #vobj <- Venn(venn.data)
    #par(bty="n")
    
    ven.fig.pep <- paste(outdir,"/all_sample-pepvenn.png",sep="",collapse="")
    png(ven.fig.pep,width=800,height=800,res=200)
    #plot(vobj, doWeights = FALSE)
    
    par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
    grob.list <- venn.diagram(x=venn.data,filename=NULL,cex=1.1,
                              cat.cex=1.4,reverse=TRUE,resolution=200,
                              height=700,width=700,compression="none")
    grid.draw(grob.list)
    dev.off()
    
    ven.fig.pep <- paste(basename(outdir),basename(ven.fig.pep),sep="/")
    
    ## 2. protein level
    venn2.data <- dlply(res$res_sample_level,
                        .variables=c(.INPUT.COLS["SAMPLE"]),function(y){
      pro <- read.delim(y[,.RES.FRACTION.COLS["PRO_SUMMARY"]],
                        colClasses=.ProSummaryColClass)
      pro <- unique(pro[,'Accession'])      
    })
    
    
    #vobj2 <- Venn(venn2.data)
    #par(bty="n")
    ven.fig.pro <- paste(outdir,"/all_sample-provenn.png",sep="",collapse="")
    png(ven.fig.pro,width=800,height=800,res=200)
    #plot(vobj2, doWeights = FALSE)
    par(mai=c(1.5,1.5,1.5,1.5),mar=c(1.5,1.5,1.5,1.5))
    grob.list <- venn.diagram(x=venn2.data, filename=NULL,
                              cex=1.1, cat.cex=1.4, reverse=TRUE,
                              resolution=200, height=700, width=700,
                              compression="none")
    grid.draw(grob.list)
    dev.off()
    
    ven.fig.pro <- paste(basename(outdir),basename(ven.fig.pro),sep="/")
    result <- data.frame(pepfig=ven.fig.pep,profig=ven.fig.pro)
  } else {
    result <- data.frame(pepfig=NA,profig=NA)
  }
  return(result)
}



plotFractionAccumulation=function(res){
  if(res$input_parameter$maxFraction >=2){
    ## For each replicate experiment
    fasta <- res$input_parameter$fasta
    outdir <- res$input_parameter$report.dir ## report figure outdir
    txtdir <- res$input_parameter$result.dir ## txt data file outdir
    ##x always is a data.frame
    tmp_para<-ddply(res$res_fraction_level,
                    #.variables=c(.INPUT.COLS["SAMPLE"],.INPUT.COLS["BIOREP"],
                    #.INPUT.COLS["TECHREP"]),
                    ## must use the following line code
                    .variables=as.vector(c(.INPUT.COLS["SAMPLE"],
                                           .INPUT.COLS["BIOREP"],
                                           .INPUT.COLS["TECHREP"])),
                    function(x){
                      ## sort the fraction from min to max
                      x <- x[order(x[,.INPUT.COLS['FRACTION']]),]
                      tmp_result <- data.frame()
                      for(i in 1:dim(x)[1]){
                        pepFile <- paste(x[1:i ,
                                           .RES.FRACTION.COLS['PEP_SUMMARY']], 
                                         sep="", collapse=";")
                        prjname=paste(paste(x[i,c(.INPUT.COLS["SAMPLE"],
                                                  .INPUT.COLS["BIOREP"],
                                                  .INPUT.COLS["TECHREP"])],
                                            sep="",collapse="_"),"_m",i,
                                      collapse="_",sep="")
                        cat("process file: ",prjname,"\n")
                        logfile = paste(txtdir,"/",prjname,"_m",i,"-log.txt",
                                        collapse="",sep="");
                        tmp_result=rbind(tmp_result,
                                         combineRun(pepFiles=pepFile,
                                                    fasta=fasta,
                                                    outPathFile=logfile,
                                                    outdir=txtdir,
                                                    prefix=prjname))
                      }
                      ## spectrum
                      fig.prefix <- paste(x[i,c(.INPUT.COLS["SAMPLE"],
                                                .INPUT.COLS["BIOREP"],
                                                .INPUT.COLS["TECHREP"])],
                                          sep="",collapse="_")
                      fig.spectra<-paste(outdir,"/",fig.prefix,"-m_spectra.png",
                                           sep="",collapse="")
                      png(fig.spectra,width=800,height=800,res=200)
                      #par(mar=c(3,3,1,1),mgp=c(1.6,0.6,0))
                      spc.df <- data.frame(x=x[,.INPUT.COLS['FRACTION']],
                                  y=tmp_result[,.RES.FRACTION.COLS['SPECTRUM']])
                      gg.obj <- ggplot(spc.df,aes(x=as.factor(x),y=y)) +
                        geom_line()+
                        geom_point()+
                        xlab("Fraction")+
                        ylab("Spectrum")
                      if(dim(x)[1] > 6){
                        gg.obj <- gg.obj + 
                          theme(axis.text.x=element_text(angle=90,vjust=0.5 ))
                        
                      }
                      print(gg.obj)
                      dev.off()
                      fig.spectra<-paste(basename(outdir),basename(fig.spectra),
                                           sep="/")

                      fig.peptide <- paste(outdir,"/",
                                           fig.prefix,"-m_peptide.png",
                                           sep="",collapse="")
                      png(fig.peptide,width=800,height=800,res=200)
                      pep.df <- data.frame(x=x[,.INPUT.COLS['FRACTION']],
                                  y=tmp_result[,.RES.FRACTION.COLS['PEPTIDE']])
                      gg.obj <- ggplot(pep.df,aes(x=as.factor(x),y=y)) +
                        geom_line()+
                        geom_point()+
                        xlab("Fraction")+
                        ylab("Peptide")
                      if(dim(x)[1] > 6){
                        gg.obj <- gg.obj + 
                          theme(axis.text.x=element_text(angle=90,vjust=0.5 ))
                        
                      }
                      print(gg.obj)
                      dev.off()
                      fig.peptide<-paste(basename(outdir),basename(fig.peptide),
                                           sep="/")
 
                      fig.protein <- paste(outdir,"/",
                                           fig.prefix,"-m_protein.png",
                                           sep="",collapse="")
                      png(fig.protein,width=800,height=800,res=200)
                      pro.df <- data.frame(x=x[,.INPUT.COLS['FRACTION']],
                                  y=tmp_result[,.RES.FRACTION.COLS['PROTEIN']])
                      gg.obj <- ggplot(pro.df,aes(x=as.factor(x),y=y)) +
                        geom_line()+
                        geom_point()+
                        xlab("Fraction")+
                        ylab("Protein")
                      if(dim(x)[1] > 6){
                        gg.obj <- gg.obj + 
                          theme(axis.text.x=element_text(angle=90,vjust=0.5 ))
                        
                      }
                      print(gg.obj)
                      dev.off()
                      fig.protein<-paste(basename(outdir),basename(fig.protein),
                                           sep="/")
                      y <- data.frame(fig.spectra=fig.spectra,
                                      fig.peptide=fig.peptide,
                                      fig.protein=fig.protein)
                      y
                      
                    })
    
    res$res_cumFraction <- tmp_para
    return(res$res_cumFraction)
  }else{
    stop("You don't have fraction design!")
  }
  
}


msQC.barplot=function(dat,xlab,ylab,fig=NULL){
  
  png(fig,width=800,height=800,res=200)
  dat$label <- sprintf("%.2f%%",100*dat$y/sum(dat$y))   
  gg.obj <- ggplot(data=dat,aes(x=as.factor(x),y=y,ymax=1.1*max(y),
                                     fill=as.factor(x))) +
    geom_bar(stat="identity",width=0.5)+
    xlab(xlab)+
    ylab(ylab)+
    scale_fill_discrete(guide=FALSE)+
    geom_text(aes(label=label),vjust=-0.2,size=3.5)
  print(gg.obj)
  dev.off()
}

