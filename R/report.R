
#############################
### report generator

##' @title HTML format report generator
##' @description HTML format report generator
##' @param res An object returned by \code{\link{msQCpipe}} function
##' @return null
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' zpqc <- system.file("extdata/qc.zip", package = "proteoQC")
##' unzip(zpqc)
##' qcres <- loadmsQCres("./qc")
##' html <- reportHTML(qcres)
reportHTML <- function(res) {

    if (is.character(res)) {
        out.report <- file.path(res, "qc_report.html")
    } else {
        out.report <- file.path(res$input_parameter$outdir, "qc_report.html")
    }
    
    if (file.exists(out.report))
        return(out.report)
    
    
    
    assign(x=".PepSummaryColClass",value=.PepSummaryColClass,envir=.GlobalEnv)
    assign(x=".ProSummaryColClass",value=.ProSummaryColClass,envir=.GlobalEnv)
    if (is.character(res) && file.exists(file.path(res, "msQC.rds")))
        res <- loadmsQCres(res)

    dir.create(res$input_parameter$report.dir,recursive=TRUE,showWarnings=FALSE)
    report<-newCustomReport("MSMS Report")    
    report<-setReportTitle(report,
                           "Quality Control Analysis of MS/MS Proteomics Data")

    
    ##The first section
    s1 <-newSection("Introduction")
    s1 <-addTo(s1,newParagraph(
        "Data quality is a concern in proteomics experiments. 
        In this report, we assess the intrinsic features of the data set 
        at multiple levels."))
   
    ##The second section: method and experiment design
    ##Firstly, there is a table to show the experiment design, 
    ##including the information of sample, replicate and fraction. 
    s2<-newSection("Methods and Data")
    ##Introduce the experiment design and data analysis
    s2_sub1 <- newSubSection(asStrong("Summary of Data Set"))
    
    design  <- ddply(res$res_fraction_level,
                     .variables=c(.INPUT.COLS['SAMPLE'],
                                  .INPUT.COLS['BIOREP'],
                                  .INPUT.COLS['TECHREP']),
                     function(x){length(x[,.INPUT.COLS['FRACTION']])})
    colnames(design) <- c(.INPUT.COLS['SAMPLE'],
                          .INPUT.COLS['BIOREP'],
                          .INPUT.COLS['TECHREP'],"No. of fraction")
    s2_sub1<-addTo(s2_sub1,newTable(design,"The Design of Experiment"))

    ##Introduce the method of data analysis and parameters.
    s2_sub2 <- newSubSection(asStrong("Data Analysis"))
    
    ##The parameters of database search 
    searchParameter <- data.frame(Item="Search engine",Value="X!Tandem")
    enzyme<-getEnzyme()$name[res$input_parameter$enzyme] ##
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Enzyme",Value=enzyme))
    fixmods <- "None"
    if(sum(res$input_parameter$fixmod %in% getMods()$index)>=1){
        fixmods <- paste(getMods()$name[res$input_parameter$fixmod],
                         collapse=";",sep="")
    }
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Fixed modifications",
                                        Value=fixmods))
    varmods <- "None"
    if(sum(res$input_parameter$varmod %in% getMods()$index)>=1){
        varmods <- paste(getMods()$name[res$input_parameter$varmod],
                         collapse=";",sep="")
    }
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Variable modifications",
                                        Value=varmods))  
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Database",
                                    Value=basename(res$input_parameter$fasta)))
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Missed cleavages",
                                  Value=as.character(res$input_parameter$miss)))
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Precursor mass error",
                                        Value=paste(res$input_parameter$tol,
                                                    res$input_parameter$tolu,
                                                    sep=" ",collapse=" ")))
    searchParameter <- rbind(searchParameter,
                             data.frame(Item="Fragment mass error",
                                        Value=paste(res$input_parameter$itol,
                                                    res$input_parameter$itolu,
                                                    sep=" ",collapse=" ")))
    searchParameterTable <- newTable(searchParameter,
                                     "The database search parameters")
    s2_sub2 <- addTo(s2_sub2,newParagraph("X!Tandem was used for analyzing the data.",
                    " Parameters used in the X!Tandem search are shown in ",
                    asReference(element=searchParameterTable),
                      ". Protein identifications were inferred from peptide identifications," ,
                      " and each identified protein had at least one associated unique peptide sequence identified at q-value equal or less than 0.01 (equivalent to a 1% FDR).", 
                      " The Occam&apos;s razor approach (Nesvizhskii, et al., 2003) was applied to deal with degenerate peptides by finding a minimum subset of proteins that covered all of the identified peptides."))
    s2_sub2<-addTo(s2_sub2,searchParameterTable)
    s2 <- addTo(s2, s2_sub1,s2_sub2)
    
  
    ##The third section: basic result 
    s3<-newSection("Results")
    ##The first subsection
    s3_sub1 <- newSubSection(asStrong("Summary of protein identification"))
    s3_sub1 <- addTo(s3_sub1,newParagraph("This part contains the basic statistics of MS/MS data."))
    ##The result of each fraction
    dview <- res$res_fraction_level[,c(.INPUT.COLS["SAMPLE"],
                                       .INPUT.COLS["BIOREP"],
                                       .INPUT.COLS["TECHREP"],
                                       .INPUT.COLS["FRACTION"],
                                       .RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                       .RES.FRACTION.COLS["SPECTRUM"],
                                       .RES.FRACTION.COLS["PEPTIDE"],
                                       .RES.FRACTION.COLS["PROTEIN"])]
    psmRatio <- sprintf("%.3f%%",100.0*dview$spectrum/dview$spectrum_total)
    dview$spectrum <- paste(dview$spectrum,"(",psmRatio,")",sep="")
    s3_sub1<-addTo(s3_sub1,newTable(dview,"Protein identification results for each fraction"))
    ##The result of each replicate experiment
    if(res$input_parameter$maxFraction>=2){
        dview <- res$res_techRep_level[,c(.INPUT.COLS["SAMPLE"],
                                          .INPUT.COLS["BIOREP"],
                                          .INPUT.COLS["TECHREP"],
                                          .RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                          .RES.FRACTION.COLS["SPECTRUM"],
                                          .RES.FRACTION.COLS["PEPTIDE"],
                                          .RES.FRACTION.COLS["PROTEIN"])]
        psmRatio <- sprintf("%.3f%%",100.0*dview$spectrum/dview$spectrum_total)
        dview$spectrum <- paste(dview$spectrum,"(",psmRatio,")",sep="")
        s3_sub1<-addTo(s3_sub1,newTable(dview,
              "Protein identification result for each technical replicate"))
    }
    ##The result of each biological replicate experiment
    if(res$input_parameter$maxTechRep>=2){
        
        dview <- res$res_bioRep_level[,c(.INPUT.COLS["SAMPLE"],
                                         .INPUT.COLS["BIOREP"],
                                         .RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                         .RES.FRACTION.COLS["SPECTRUM"],
                                         .RES.FRACTION.COLS["PEPTIDE"],
                                         .RES.FRACTION.COLS["PROTEIN"])]
        psmRatio <- sprintf("%.3f%%",100.0*dview$spectrum/dview$spectrum_total)
        dview$spectrum <- paste(dview$spectrum,"(",psmRatio,")",sep="")
        s3_sub1<-addTo(s3_sub1,newTable(dview,
              "Protein identification result for each biological replicate"))
    }

    ##The result of each sample
    if(res$input_parameter$maxBioRep>=2){
  
        dview <- res$res_bioRep_level[,c(.INPUT.COLS["SAMPLE"],
                                         .RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                         .RES.FRACTION.COLS["SPECTRUM"],
                                         .RES.FRACTION.COLS["PEPTIDE"],
                                         .RES.FRACTION.COLS["PROTEIN"])]
        psmRatio <- sprintf("%.3f%%",100.0*dview$spectrum/dview$spectrum_total)
        dview$spectrum <- paste(dview$spectrum,"(",psmRatio,")",sep="")
        s3_sub1<-addTo(s3_sub1,newTable(dview,
              "Protein identification result for each sample"))
    }
    s3 <- addTo(s3,s3_sub1)
    
    summary.chart.subsection <- newSubSection(asStrong("Summary charts"))
    summary.chart.subsection <- addTo(summary.chart.subsection,
                                      addSummaryChart(res))
    s3 <- addTo(s3,summary.chart.subsection)
    
    ## contaminants stat
    
    cnt.stat.subsetion <- newSubSection(asStrong("Contaminants stat"))
    cnt.stat.subsetion <- addTo(cnt.stat.subsetion,
        newParagraph("The common Repository of Adventitious Proteins, ",
                    "cRAP (pronounced \"cee-RAP\"), ",
                    "is an attempt to create a list of proteins commonly found 
                    in proteomics experiments that are present either by 
                    accident or through unavoidable contamination of protein 
                    samples. ",
                    "The types of proteins included fall into three general 
                    classes:"),
        newList(isNumbered=TRUE,  
                newParagraph("common laboratory proteins;"),
                newParagraph(
        "proteins added by accident through dust or physical contact; and"),
                newParagraph(
        "proteins used as molecular weight or mass spectrometry quantitation 
        standards.")),
                newParagraph("We added the ",
                             asLink("http://www.thegpm.org/crap/"," cRAP"),
                             " database in database searching."))
    
    
    if (!is.null(cntStat(res))) {
        cnt.stat.subsetion <- addTo(cnt.stat.subsetion,
                                  newTable(table=cntStat(res),
                                    "Identification of contaminant proteins"))
        s3 <- addTo(s3,cnt.stat.subsetion)
    } else {
        cnt.stat.subsetion <- addTo(cnt.stat.subsetion,newParagraph(
                              "There is none contaminant protein identified!"))
        s3 <- addTo(s3,cnt.stat.subsetion)
    }
    
    
    ##Reproducibility evaluation
    s3_sub2 <- newSubSection(asStrong("Reproducibility"))
    
    
    s3_sub2_sub1 <- 
      newSubSubSection("Reproducibility of result for each fraction")
    
    tmp_sub <- newParagraph(
      asStrong("Reproducibility of total spectra for each fraction"))

    specTotal<-plotFractionIDResult(res,level=1)

    specTotal_errorbar<-plotSampleIDResultErrorBar(res,level=1)

    tmp_sub <- addTo(tmp_sub,
          newFigure(specTotal,"The distribution of the total spectra number"))
    tmp_sub <- addTo(tmp_sub,
          newFigure(specTotal_errorbar,
              "Error bar plot of the total spectra number for each fraction"))
    s3_sub2_sub1 <- addTo(s3_sub2_sub1 ,tmp_sub)

    tmp_sub <- newParagraph(
      asStrong("Reproducibility of identified spectra for each fraction"))

    specTotal<-plotFractionIDResult(res,level=2)

    specTotal_errorbar<-plotSampleIDResultErrorBar(res,level=2)

    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal,
                          "The distribution of the identified spectra number"))
    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal_errorbar,
                  "Error bar plot of the identified spectra for each fraction"))

    s3_sub2_sub1 <- addTo(s3_sub2_sub1,tmp_sub)
    tmp_sub <- newParagraph(
      asStrong("Reproducibility of identified peptides for each fraction"))
 
    specTotal<-plotFractionIDResult(res,level=3)
  
    specTotal_errorbar<-plotSampleIDResultErrorBar(res,level=3)
    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal,
                            "The distribution of identified peptides number"))
    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal_errorbar,
                "Error bar plot of the identified peptides for each fraction"))

    s3_sub2_sub1 <- addTo(s3_sub2_sub1,tmp_sub)
    
    tmp_sub <- newParagraph(
      asStrong("Reproducibility of identified proteins for each fraction"))
  
    specTotal<-plotFractionIDResult(res,level=4)
   
    specTotal_errorbar<-plotSampleIDResultErrorBar(res,level=4)
    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal,
                        "The distribution of the identified proteins number"))
    tmp_sub<-addTo(tmp_sub,
                   newFigure(specTotal_errorbar,
                "Error bar plot of the identified proteins for each fraction"))
    
    
    s3_sub2_sub1 <- addTo(s3_sub2_sub1,tmp_sub)
    s3_sub2 <- addTo(s3_sub2,s3_sub2_sub1)

    ##The eeproducibility of technical replicates
    ##peptide/protein
    if(res$input_parameter$maxTechRep>=2){
        ## cat("12\n")
        s3_sub2_sub2 <- 
          newSubSubSection("Reproducibility of result for technical replicate")
        venn.fig <- plotTechRepVenn(res)
        
        ## venn.fig
        for(i in 1:dim(venn.fig)[1]){
            x <- venn.fig[i,]
            ## paste(x[1,"SAMPLE"],x[1,"BIOREP"],sep=" ",collapse=)
            
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
              newFigure(x[1,"pepfig"],
                  "Technical replicate of peptides' overlap between sample ",
                  x[1,"SAMPLE"]," biological replicate ",x[1,"BIOREP"],"."))
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
              newFigure(x[1,"profig"],
                  "Technical replicate of proteins' overlap between sample ",
                  x[1,"SAMPLE"]," biological replicate ",x[1,"BIOREP"],"."))
            
        }
        s3_sub2 <- addTo(s3_sub2,s3_sub2_sub2)
        
    }
  
    ##Reproducibility of result for biological replicate
    if(res$input_parameter$maxBioRep>=2){
        ## cat("13\n")
        s3_sub2_sub2 <- 
          newSubSubSection("Reproducibility of result for biological replicate")
        venn.fig <- plotBioRepVenn(res)
        ##venn.fig
        for(i in 1:dim(venn.fig)[1]){
            x <- venn.fig[i,]
            ## paste(x[1,"SAMPLE"],x[1,"BIOREP"],sep=" ",collapse=)
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
                newFigure(x[1,"pepfig"],
                    "Biological replicate of peptides' overlap between sample ",
                    x[1,"SAMPLE"],"."))
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
                newFigure(x[1,"profig"],
                    "Biological replicate of proteins' overlap between sample ",
                    x[1,"SAMPLE"],"."))
            
        }
        s3_sub2 <- addTo(s3_sub2,s3_sub2_sub2)
        ## s3<-addTo(s3,s3_sub2)
    }
    
    ##the overlap of samples
    if(res$input_parameter$maxSample>=2){
        ## cat("14\n")
        s3_sub2_sub2 <- 
          newSubSubSection("Reproducibility of result for different samples")
        venn.fig <- plotSampleVenn(res)
        ## venn.fig
        for(i in 1:dim(venn.fig)[1]){
            x <- venn.fig[i,]
            ## paste(x[1,"SAMPLE"],x[1,"BIOREP"],sep=" ",collapse=)
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
                            newFigure(x[1,"pepfig"],
                              "Peptides' overlap between different samples."))
            s3_sub2_sub2<-addTo(s3_sub2_sub2,
                            newFigure(x[1,"profig"],
                              "Proteins' overlap between different samples."))
            
        }
        s3_sub2 <- addTo(s3_sub2,s3_sub2_sub2)
        ## s3<-addTo(s3,s3_sub2)
    }
    s3<-addTo(s3,s3_sub2)
    
    
    ##Mass accuracy
    s3_sub3 <- newSubSection(asStrong("Mass accuracy"))
    
    
    ms2.fig <- plotMS2Error(res)
    
    ms1ppm.fig <- plotMS1Error(res,plot.class="ppm")
    
    ms1da.fig <- plotMS1Error(res,plot.class="da")
    
    s3_sub3 <- addTo(s3_sub3,
                     newFigure(ms2.fig,
                               "The mass error (Da) of the fragment ions"))
    s3_sub3 <- addTo(s3_sub3,
                     newFigure(ms1ppm.fig,
                               "The mass error (ppm) of the precusor"))
    s3_sub3 <- addTo(s3_sub3,
                     newFigure(ms1da.fig,
                               "The mass error (Da) of the precusor"))
    s3 <- addTo(s3,s3_sub3)
    ##Separation efficiency
    ##
    s3_sub4 <- newSubSection(asStrong("Separation efficiency"))
    
    if(res$input_parameter$maxFraction>=2){
        fig.cum <- plotFractionAccumulation(res)
        s3_sub4_sub1 <- newSubSubSection("Fractionation cumulative effect")
        for(i in 1:dim(fig.cum)[1]){
            x <- fig.cum[i,]
            fig.name.prefix <- paste("Sample: ",
                                     x[1,.INPUT.COLS["SAMPLE"]],
                                     ", biological replicate: ",
                                     x[1,.INPUT.COLS["BIOREP"]],
                                     ", technical replicate: ",
                                     x[1,.INPUT.COLS["TECHREP"]],
                                     sep=" ",collapse=" ")
            s3_sub4_sub1<-addTo(s3_sub4_sub1,
                                newFigure(x[1,"fig.spectra"],
                                  fig.name.prefix,
                                  ", identified spectrum accumulative effect."))
            s3_sub4_sub1<-addTo(s3_sub4_sub1,
                                newFigure(x[1,"fig.peptide"],
                                  fig.name.prefix,
                                  ", identified peptides accumulative effect."))
            s3_sub4_sub1<-addTo(s3_sub4_sub1,
                                newFigure(x[1,"fig.protein"],
                                  fig.name.prefix,
                                  ", identified proteins accumulative effect."))
            
        }
        s3_sub4<-addTo(s3_sub4,s3_sub4_sub1)
        
    }
    s3<-addTo(s3,s3_sub4)

    ## Generate identification-independent metrics
    if(!is.null(res$ms12metrics)){
      outfig <- plotMS12(res$ms12metrics,res$input_parameter$report.dir)  
      IDfree_sub <- newSubSection(
                        asStrong("Identification-independent quality metrics"))
      
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1tic),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "TIC of MS1 distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1peakscount),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "MS1 peak count distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1ioncount),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "MS1 ion count distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms2peaksdensity),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                              "MS2 peak density distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms2boxplot),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "MS2 peak boxplot"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1countdot),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "MS1 peak count distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1counterrorbar),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,
                                                "MS1 peak count distribution"))
      fig <- paste(basename(res$input_parameter$report.dir),
                   basename(outfig$ms1boxplot),
                   sep="/")
      IDfree_sub <- addTo(IDfree_sub, newFigure(fig,                                          
                                                "MS1 peak boxplot"))
      s3<-addTo(s3,IDfree_sub)
    }
    
    report<-addTo(report,s1,s2,s3)
    fn <- file.path(res$input_parameter$outdir,"qc_report")
    writeReport(report, filename = fn)
    fn <- paste0(fn, ".html")
    message("QC report written to ", fn)
    invisible(fn)
}


##' @title Add PRIDE summary charts
##' @description Add PRIDE summary charts in the technical replicate level
##' @param res An object returned by msQCpipe function
##' @return NULL
addSummaryChart <- function(res){  
  x <- NULL
  outdir <- res$input_parameter$report.dir
  par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1))
  if(res$input_parameter$maxFraction==1){
    # no fraction
    x <- res$res_fraction_level
  }else{
    # have fraction
    x <- res$res_techRep_level
  }
  
  out.report <- 
    newParagraph("Summary plot for each technical replicate experiment.")
  report.prefix.dir <- basename(outdir)
  
  for(i in 1:dim(x)[1]){
    prefix <- paste(x[i,c(.INPUT.COLS['SAMPLE'],
                          .INPUT.COLS['BIOREP'],
                          .INPUT.COLS['TECHREP'])],
                    sep="",collapse="_")
    html.fig.prefix <- paste("sample: ",x[i,.INPUT.COLS['SAMPLE']],
                             ", biological replicate: ",
                             x[i,.INPUT.COLS['BIOREP']],
                             ", technical replicate: ",
                             x[i,.INPUT.COLS['TECHREP']],sep="",collapse="")
    ## html report
    html.report <- newParagraph("Summary charts for sample: ",
                                x[i,.INPUT.COLS['SAMPLE']],
                                ", biological: ",
                                x[i,.INPUT.COLS['BIOREP']],
                                ", technical: ",
                                x[i,.INPUT.COLS['TECHREP']],
                                ".");
    
    # read peptide/protein summary
    peps <- read.delim(x[i,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                       stringsAsFactors=FALSE,colClasses=.PepSummaryColClass)
    pros <- read.delim(x[i,.RES.FRACTION.COLS["PRO_SUMMARY"]],
                       stringsAsFactors=FALSE,colClasses=.ProSummaryColClass)
    # Missed Tryptic Cleavages, xtandem output bug, 
    # there are many large miss peptide
    tmp.miss <- peps$miss
    #tmp.miss[peps$miss>res$input_parameter$miss] <- res$input_parameter$miss
    tmp.miss[peps$miss>res$input_parameter$miss] <- Inf
  
    fig.name <- paste(outdir,"/",prefix,"-miss.png",sep="")

    miss.df  <- as.data.frame(table(tmp.miss))
    names(miss.df) <- c("x","y")
    msQC.barplot(miss.df,xlab="Missed cleavages",
                 ylab="Number of peptides",fig=fig.name)
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                            "Missed cleavages chart for ",html.fig.prefix,"."))
    
    ## precursor ion charge
    fig.name <- paste(outdir,"/",prefix,"-charge.png",sep="")
    charge.data <- as.data.frame(table(peps$charge))
    names(charge.data) <- c("x","y")
    msQC.barplot(charge.data,xlab="Peptide charge",
                 ylab="Number of peptides",fig=fig.name)
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                      newFigure(file=fig.name, 
                        "Precursor ion charge chart for ",html.fig.prefix,"."))
    
    ## peptide length
    pep.length <- nchar(peps$peptide)
    fig.name <- paste(outdir,"/",prefix,"-peplength.png",sep="")
    #png(fig.name,width=500,height=500)
    png(fig.name,width=1000,height=800,res=200)
#   par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1),cex.lab=1.2,cex.axis=1.1,font.lab=2)
#   par(yaxs = "i")    #plot(density(mass.data),main="",xlab="Mass (Daltons)")
#   tmp.h <- hist(pep.length,nclass=50,plot=FALSE)
#   hist(pep.length,nclass=50,main="",
#        xlab="Peptide length",ylim=c(0,1.1*max(tmp.h$counts)))
#   box()
    p.dataframe <- data.frame(peplength=pep.length)
    gg.obj <- ggplot(data=p.dataframe,mapping=aes(x=peplength))
    #gg.obj <- gg.obj + geom_histogram(binwidth=1,fill="white",colour="purple") 
    gg.obj <- gg.obj + xlab("Peptide length") +
              geom_histogram(aes(fill = ..count..),colour="white")+
              scale_fill_gradient("Count", low = "green", high = "red")
              #geom_rug(aes(y=peplength),sides="b",position="jitter",alpha=0.4,colour="gray");
    print(gg.obj)
    dev.off()
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                              "Peptide length chart for ",html.fig.prefix,"."))
    
    ## peaks per ms/ms spectrum
    ## TODO, currently, the ms number extracted from the X!Tandem 
    ## xml result file is not correct
    
    ## precursor ion error delta
    # Maybe there are some precursor ion mass delta (da) 
    # approximately equal +-1 because of isotope error
    delta.da  <- peps$delta_da
    delta.ppm <- as.numeric(peps$delta_ppm)
    if(res$input_parameter$tolu == "ppm"){
      
      delta.da  <- delta.da[abs(delta.ppm) <= res$input_parameter$tol]
      delta.ppm <- delta.ppm[abs(delta.ppm) <= res$input_parameter$tol]
      
    }else{
      ## dalton
      delta.da  <- delta.da[abs(delta.da) <= res$input_parameter$tol]
      delta.ppm <- delta.ppm[abs(delta.da) <= res$input_parameter$tol]
      
    }
    
    # Da plot
    fig.name <- paste(outdir,"/",prefix,"-ms1delta_da.png",sep="")
    #png(fig.name,width=500,height=500)
    png(fig.name,width=1000,height=800,res=200)
    delta.da.df <- data.frame(delta=delta.da)
    r3 <- quantile(delta.da,probs=c(0.05,0.50,0.95))
    #save(delta.da.df,file="ok.rda")
    r3.df <- data.frame(r3=r3)
    gg.obj <- ggplot(data=delta.da.df,mapping=aes(x=delta))+
                geom_histogram(aes(fill = ..count..),colour="white")+
                scale_fill_gradient("Count", low = "green", high = "red")+
                #geom_histogram(fill="white",colour="purple")+
                xlab("Precursor mass delta/Da")+
                geom_vline(data=r3.df,mapping=aes(xintercept=r3),
                           colour="blue",size=1,alpha=0.5) +
                annotate(geom="text",x=Inf,y=Inf,vjust=1,hjust=1,size=3,
                         label=paste("Quantiles:",sprintf("5%%=%.4f" ,r3[1]),
                                     sprintf("50%%=%.4f",r3[2]),
                                     sprintf("95%%=%.4f",r3[3]),
                                     sprintf("[5%%,95%%]=%.4f",r3[3]-r3[1]),
                                     sep="\n",collapse=""))
    print(gg.obj)
    dev.off()
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report<-addTo(html.report,
                  newFigure(file=fig.name, 
                    "Precursor mass delta (Da) chart for ",html.fig.prefix,"."))
    
    # ppm plot
    fig.name <- paste(outdir,"/",prefix,"-ms1delta_ppm.png",sep="")
    png(fig.name,width=1000,height=800,res=200)
    delta.da.df <- data.frame(delta=delta.ppm)
    r3 <- quantile(delta.ppm,probs=c(0.05,0.50,0.95))
    #save(delta.da.df,file="ok.rda")
    r3.df <- data.frame(r3=r3)
    gg.obj <- ggplot(data=delta.da.df,mapping=aes(x=delta))+
      geom_histogram(aes(fill = ..count..),colour="white")+
      scale_fill_gradient("Count", low = "green", high = "red")+
      #geom_histogram(fill="white",colour="purple")+
      xlab("Precursor mass delta/ppm")+
      geom_vline(data=r3.df,mapping=aes(xintercept=r3),
                 colour="blue",size=1,alpha=0.5) +
      annotate(geom="text",x=Inf,y=Inf,vjust=1,hjust=1,size=3,
               label=paste("Quantiles:",sprintf("5%%=%.4f" ,r3[1]),
                           sprintf("50%%=%.4f",r3[2]),
                           sprintf("95%%=%.4f",r3[3]),
                           sprintf("[5%%,95%%]=%.4f",r3[3]-r3[1]),
                           sep="\n",collapse=""))
    print(gg.obj)
    dev.off()
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                            "Precursor mass delta (ppm) chart for ",
                            html.fig.prefix,"."))

    
    ## fragment ion mass delta


    ms2delta.data <- as.numeric(unlist(strsplit(x=as.character(peps$ms2delta),
                                                split=";")))
    fig.name <- paste(outdir,"/",prefix,"-ms2delta_da.png",sep="")
    png(fig.name,width=1000,height=800,res=200)
    delta.da.df <- data.frame(delta=ms2delta.data)
    r3 <- quantile(ms2delta.data,probs=c(0.05,0.50,0.95))
    #save(delta.da.df,file="ok.rda")
    r3.df <- data.frame(r3=r3)
    gg.obj <- ggplot(data=delta.da.df,mapping=aes(x=delta))+
      geom_histogram(aes(fill = ..count..),colour="white")+
      scale_fill_gradient("Count", low = "green", high = "red")+
      #geom_histogram(fill="white",colour="purple")+
      xlab("Fragment ion mass delta/Da")+
      geom_vline(data=r3.df,mapping=aes(xintercept=r3),
                 colour="blue",size=1,alpha=0.5) +
      annotate(geom="text",x=Inf,y=Inf,vjust=1,hjust=1,size=3,
               label=paste("Quantiles:",sprintf("5%%=%.4f" ,r3[1]),
                           sprintf("50%%=%.4f",r3[2]),
                           sprintf("95%%=%.4f",r3[3]),
                           sprintf("[5%%,95%%]=%.4f",r3[3]-r3[1]),
                           sep="\n",collapse=""))
    print(gg.obj)
    dev.off()
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                                   "Fragment ion mass delta (Da) chart for ",
                                   html.fig.prefix,"."))
    
    
    ## retention time plot
    if(is.numeric(peps$rt) && peps$rt[1]>=0){
      rt.data <- peps$rt/60
      fig.name <- paste(outdir,"/",prefix,"-rt.png",sep="")
      png(fig.name,width=500,height=500)
      par(mgp=c(1.6,0.6,0),mar=c(3,3,1,1),cex.lab=1.2,cex.axis=1.1,font.lab=2)
      par(yaxs = "i")    #plot(density(mass.data),main="",xlab="Mass (Daltons)")
      tmp.h <- hist(rt.data,nclass=50,plot=FALSE)
      tmp.h <- hist(rt.data,nclass=50,plot=FALSE)
      hist(rt.data,nclass=50,main="",
           xlab="Retention time",
           ylim=c(0,1.1*max(tmp.h$counts)))
      lines(density(rt.data))
      box()
      dev.off()
      fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
      html.report <- addTo(html.report,
                           newFigure(file=fig.name, 
                                     "Retention time chart for ",
                                     html.fig.prefix,"."))
    }
    
    ## protein level plot
    
    ## unique spectrum per protein
    #save(pros,file="pros.rda")
    pro.ms2 <- pros$NumOfUniqSpectra
    ms2.data <- table(cut(pro.ms2,
                          breaks=c(0,1,2,3,5,10,Inf),
                          labels=c("1","2","3","(3,5]","(5,10]",">10")))
    fig.name <- paste(outdir,"/",prefix,"-uniqueSpectrumPerProtein.png",sep="")
    bar.labels <- sprintf("%.2f%%",100*ms2.data/sum(ms2.data))   
    spc.dat <- as.data.frame(ms2.data)
    names(spc.dat) <- c("x","y")
    msQC.barplot(spc.dat,xlab="Number of spectrum",
                 ylab="Number of proteins",fig=fig.name)
    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                                   "Unique spectrum per protein chart for ",
                                   html.fig.prefix,"."))
    
    ## unique peptides per protein
    pro.npep <- pros$NumOfUniqPeps
    npep.data <- table(cut(pro.npep,breaks=c(0,1,2,3,5,10,Inf),
                           labels=c("1","2","3","(3,5]","(5,10]",">10")))
    fig.name <- paste(outdir,"/",prefix,"-uniquePeptidePerProtein.png",sep="")
    npep.dat <- as.data.frame(npep.data)
    names(npep.dat) <- c("x","y")
    msQC.barplot(npep.dat,xlab="Number of peptides",
                 ylab="Number of proteins",fig=fig.name)

    fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
    html.report <- addTo(html.report,
                         newFigure(file=fig.name, 
                                   "Unique peptide per protein chart for ",
                                   html.fig.prefix,"."))
    ## protein masses
    if(is.numeric(pros$Mass)){
      pro.mass <- pros$Mass/1000
      fig.name <- paste(outdir,"/",prefix,"-promass.png",sep="")
      png(fig.name,width=1000,height=800,res=200)
      proMass.df <- data.frame(pro.mass=pro.mass)
      gg.obj <- ggplot(data=proMass.df,mapping=aes(x=pro.mass))
      #gg.obj <- gg.obj + geom_histogram(binwidth=1,fill="white",colour="purple") 
      gg.obj <- gg.obj + xlab("Protein mass (kDa)") +
        geom_histogram(aes(fill = ..count..),colour="white")+
        scale_fill_gradient("Count", low = "green", high = "red")
      #geom_rug(aes(y=peplength),sides="b",position="jitter",alpha=0.4,colour="gray");
      print(gg.obj)
      
      dev.off()
      fig.name <- paste(report.prefix.dir,"/",basename(fig.name),sep="")
      html.report <- addTo(html.report,
                           newFigure(file=fig.name, 
                                     "Protein mass chart for ",
                                     html.fig.prefix,"."))
    }
    
    out.report <- addTo(out.report,html.report)
    
  }
  return(out.report)
}

##' @title contaminants stat
##' @description Common Contaminants in Proteomics Mass Spectrometry Experiments
##' @param res An object of msQCres
##' @return A data.frame will be shown in HTML report
cntStat=function(res){
  
  
  cnt.map <- read.delim(file=system.file("db/crap.txt",package="proteoQC"),
                        stringsAsFactors=FALSE)
  
  ##
  x <- NULL
  if(res$input_parameter$maxFraction ==1){
    x <- res$res_fraction_level
  }else{
    x <- res$res_techRep_level
  }
  
  y <- data.frame()
  for(i in 1:dim(x)[1]){
    pros <- read.delim(x[i,.RES.FRACTION.COLS["PRO_SUMMARY"]],
                       colClasses=.ProSummaryColClass)
    cnt.accs <- pros[grep(pattern="cnt_msQC_",x=pros$Accession),]
    if(dim(cnt.accs)[1] >=1){
      tmp <- cnt.accs[,c("Accession","NumOfUniqPeps","NumOfUniqSpectra"),
                      drop=FALSE]
      tmp$Sample <- rep(paste(x[i,c(.INPUT.COLS["SAMPLE"],
                                    .INPUT.COLS["BIOREP"],
                                    .INPUT.COLS["TECHREP"])],
                              sep="",collapse="_"),dim(tmp)[1])
      y <- rbind(y,tmp)
    }
  }


  ## It's possible that there are none contamitant protein identified, 
  ## so maybe y is none value
  if(nrow(y)>=1){
    y$Accession <- sub(pattern="cnt_msQC_",replacement="",x=y$Accession)
    names(y) <- c("Accession","Peptides","Spectrum","Sample")
    yy <- merge(y,cnt.map,by.x="Accession",by.y="Acc",sort=TRUE)
    if(dim(yy)[1]!=dim(y)[1]){
      warnings("Some contaminants protein are not in crap.txt!\n")
    }
    yy <- yy[order(yy$Sample,yy$Spectrum,decreasing=TRUE),]
  }else{
    yy<-NULL
  }
  return(yy)
}


