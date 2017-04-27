
##' @title The main function of msQC pipeline
##'
##' @description This function is designed to automate generating of 
##' target-decoy database, database searcing, post-processing and report 
##' generation.
##' @param spectralist A file contains the experiment design or a single mgf file
##' @param fasta database file, must contain decoy sequences
##' @param outdir output directory
##' @param mode identification or quantification
##' @param miss max miss clevage
##' @param enzyme enzyme
##' @param varmod Variable modifications are those which may or may
##' not be present.
##' @param fixmod Fixed modifications are applied universally, to
##' every instance of the specified residue(s) or terminus.
##' @param tol The error window on experimental peptide mass values
##' @param tolu Units can be selected from: ppm, Daltons(also da or Da).
##' @param itol Error window for MS/MS fragment ion mass values.
##' @param itolu  Units can be selected from: Daltons(also da or Da)
##' @param refine Refine search for X!Tandem, default is TRUE.
##' @param ntt Semi-tryptic, 1; fully-tryptic, 2.
##' @param threshold  FDR value for PSM
##' @param cpu  Max number of cpu used
##' @param xmx  JAVA -Xmx
##' @param ...  Additional parameters passed to
##' \code{\link{read.table}} used to read the experimental design.
##' @return A list which contains all of the information for data 
##' quality report generating
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' \dontrun{
##' library("rpx")
##' px <- PXDataset("PXD000864")
##' mgfs <- grep("mgf", pxfiles(px), value = TRUE)
##' mgfs <- grep("-0[5-6]-[1|2]", mgfs, value=TRUE)
##' mgffiles <- pxget(px, mgfs)
##' library("R.utils")
##' mgffiles <- sapply(mgffiles, gunzip)
##' ## Generate the lightweight qc report, 
##' ## trim the mgf files to 1/10 of their size.
##' trimMgf <- function(f, m = 1/10, overwrite = FALSE) {
##'   message("Reading ", f)
##'   x <- readLines(f)
##'   beg <- grep("BEGIN IONS", x)
##'   end <- grep("END IONS", x)
##'   n <- length(beg)
##'   message("Sub-setting to ", m)
##'   i <- sort(sample(n, floor(n * m)))
##'   k <- unlist(mapply(seq, from = beg[i], to = end[i]))
##'   if (overwrite) {
##'     unlink(f)
##'     message("Writing ", f)
##'     writeLines(x[k], con = f)
##'     return(f)
##'   } else {
##'     g <- sub(".mgf", "_small.mgf", f)
##'     message("Writing ", g)
##'     writeLines(x[k], con = g)
##'     return(g)
##'   }    
##' }
##' set.seed(1)
##' mgffiles <- sapply(mgffiles, trimMgf, overwrite = TRUE)
##' fas <- pxget(px, "TTE2010.zip")
##' fas <- unzip(fas)
##' design <- system.file("extdata/PXD000864-design.txt", package = "proteoQC")
##' read.table(design, header = TRUE)
##' qcres <- msQCpipe(spectralist = design,
##'                  fasta = fas, 
##'                  outdir = "./qc",
##'                  miss  = 0,
##'                  enzyme = 1, varmod = 2, fixmod = 1,
##'                  tol = 10, itol = 0.6, cpu = 2,
##'                  mode = "identification")
##' }
msQCpipe <- function(spectralist=NULL, fasta="", outdir="./", mode="",
                     miss=2, enzyme=1,
                     varmod=NULL, fixmod=NULL, ##modification
                     tol=10, tolu="ppm", itol=0.6, itolu="Daltons", ##mass error
                     threshold=0.01, cpu=0, xmx=2,refine=TRUE,ntt=1,
                     ...) {    
  ## check if the input parameters are valid!
  if (!file.exists(spectralist)) {
    stop("The design file: ", spectralist, " does not exists!")
  }
  if(tolu!="ppm"){
    tolu="Daltons"
  }
  if(itolu!="ppm"){
    itolu="Daltons"
  }
    
  ## take an mgf file as input is supported.
  if(grepl(pattern = ".mgf",x = spectralist,ignore.case = TRUE)){
    ## generate a design file in outdir.
    dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
    tmp_spectralist <- data.frame(file=spectralist,sample=1, bioRep=1, techRep=1, fraction=1)
    tmp_spectralist_file <- paste(outdir,"/sample_list.txt",sep = "",collapse = "")
    write.table(tmp_spectralist,file = tmp_spectralist_file,col.names = TRUE,
                row.names = FALSE,quote = FALSE,sep="\t")
    spectralist <- tmp_spectralist_file
  }
  
  ## save all the information we need, including the input parameters,
  ## the information needed by report generation
  res<-list(input_parameter=list())
  ## save the input parameter into the res object
  res$input_parameter$spectralist=spectralist  
  res$input_parameter$outdir=outdir
  res$input_parameter$mode=mode
  res$input_parameter$miss=miss
  res$input_parameter$enzyme=enzyme
  res$input_parameter$varmod=varmod
  res$input_parameter$fixmod=fixmod
  res$input_parameter$tol=tol
  res$input_parameter$tolu=tolu
  res$input_parameter$itol=itol
  res$input_parameter$itolu=itolu
  res$input_parameter$refine=refine
  res$input_parameter$ntt=ntt
  res$input_parameter$threshold=threshold
  res$input_parameter$date=date()
  
  
  ## Create the protein result file directory
  res$input_parameter$result.dir <- paste(res$input_parameter$outdir,"/result",
                                          sep="")
  dir.create(res$input_parameter$result.dir,recursive=TRUE,showWarnings=FALSE)
  
  ## Create the report file directory
  res$input_parameter$report.dir <- paste(res$input_parameter$outdir,"/report",
                                          sep="")
  dir.create(res$input_parameter$report.dir,recursive=TRUE,showWarnings=FALSE)
  
  ## Create the database file directory
  res$input_parameter$db.dir = paste(res$input_parameter$outdir,"/database",
                                     sep="")
  dir.create(res$input_parameter$db.dir,recursive=TRUE,showWarnings=FALSE)
  
  ## target decoy database file path
  res$input_parameter$fasta <- paste0(res$input_parameter$db.dir,
                                      "/target_decoy.fasta")
  createTargetDecoyDB(fasta,res$input_parameter$fasta)
  
  res$input_parameter$target.db=fasta
  
  ## cpu
  if (cpu == 0) {
    cpu = detectCores()
  } else {
    cpu_tmp = detectCores()
    if (cpu > cpu_tmp) {
      warning("The number of cpus is too large, set it to ", cpu, sep = " ")
      cpu = cpu_tmp
    }
    
  }
  
  res$input_parameter$cpu=cpu
  message("Using ", cpu," cpu's\n")
  
  ## read the input parameters
  para <- read.table(spectralist, header=TRUE, stringsAsFactors=FALSE, ...)
  
  res$input_parameter$maxFraction=length(unique(para$fraction))
  res$input_parameter$maxTechRep=length(unique(para$techRep))
  res$input_parameter$maxBioRep=length(unique(para$bioRep))
  res$input_parameter$maxSample=length(unique(para$sample))
  
  ## run X!Tandem for each mgf file
  result<-data.frame()
  for(i in 1:dim(para)[1]){
    mgf=para[i,.INPUT.COLS["FILE"]]
    #prjname=paste(para[i,-1],collapse="_",sep="")
    prjname <- paste(para[i,c(.INPUT.COLS["SAMPLE"], .INPUT.COLS["BIOREP"],
                              .INPUT.COLS["TECHREP"], .INPUT.COLS["FRACTION"])],
                     collapse="_", sep="")
    message("Processing file: ", prjname, "\n")
    result=rbind(result,
        runTandem(spectra=mgf,
                  outdir=res$input_parameter$result.dir,
                  outprefix=prjname,  
                  varmod=varmod, fixmod=fixmod,
                  tol=tol, tolu=tolu,
                  itol=itol, itolu=itolu,
                  fasta=res$input_parameter$fasta,
                  cpu=cpu,
                  enzyme=enzyme,
                  xmx=xmx,                                  
                  miss=miss,refine = refine,ntt = ntt))
    
  }
  
  result<-cbind(para,result)
  ## fraction level result, it's a data.frame object
  res$res_fraction_level=result
  
  ## sample, bioreplicate and techreplicate
  ## stat for different fraction of the same tech replicate
  ## We need to judge whether the ms data have more that one fraction. 
  ## If no more that one fraction for each experiment, then we need not 
  ## to do this step.
  ## Only when number of fraction > 1, then need to combine the fractions 
  ## for a replicate
  if(res$input_parameter$maxFraction>1){
    ## x is always data.frame
    #tmp_para<-ddply(res$res_fraction_level,.(sample,bioRep,techRep),
    #function(x){
    tmp_para<-ddply(res$res_fraction_level,
                #.variables=c(.INPUT.COLS["SAMPLE"],.INPUT.COLS["BIOREP"],
                #.INPUT.COLS["TECHREP"]),
                ##Must use the following code, the above code cann't work well
                .variables=as.vector(c(.INPUT.COLS["SAMPLE"],
                                       .INPUT.COLS["BIOREP"],
                                       .INPUT.COLS["TECHREP"])),
                function(x){    
                  y=data.frame(x=sum(x[,.RES.FRACTION.COLS["SPECTRUM_TOTAL"]]),
                               y=paste(x[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                                       sep="",collapse=";"))
                  colnames(y)=c(.RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                .RES.COMBINE.COLS["PEP_SUMMARY_FILES"])
                  y
                })
    tmp_result<-data.frame()
    for (i in 1:dim(tmp_para)[1]) {
        files=tmp_para[i,.RES.COMBINE.COLS["PEP_SUMMARY_FILES"]]
        prjname <- paste(tmp_para[i,c(.INPUT.COLS["SAMPLE"],
                                      .INPUT.COLS["BIOREP"],
                                      .INPUT.COLS["TECHREP"])],
                         collapse="_", sep="")
        message("processning file: ",prjname,"\n")
        logfile = paste(res$input_parameter$result.dir, "/", 
                        prjname, "-log.txt", collapse="", sep="");
        tmp_result <- rbind(tmp_result,
                            combineRun(pepFiles=files, 
                                       fasta=res$input_parameter$fasta,
                                       outPathFile=logfile, 
                                       outdir=res$input_parameter$result.dir,
                                       prefix=prjname))
    }
    
    tmp_result<-cbind(tmp_para,tmp_result)
    res$res_techRep_level=tmp_result
  } else {
    ## if no fraction , then 
    res$res_techRep_level <- res$res_fraction_level
  }
  
  
  
  ## Combine result for the biological replicate
  if(res$input_parameter$maxTechRep>1){
    ## x is always data.frame
    tmp_para<-ddply(res$res_fraction_level,
                .variables=as.vector(c(.INPUT.COLS["SAMPLE"],
                                       .INPUT.COLS["BIOREP"])),
                function(x){
                  y=data.frame(x=sum(x[,.RES.FRACTION.COLS["SPECTRUM_TOTAL"]]),
                               y=paste(x[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                                       sep="",collapse=";"))
                  colnames(y)=c(.RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                .RES.COMBINE.COLS["PEP_SUMMARY_FILES"])
                  y
                })
    tmp_result<-data.frame()
    for(i in 1:dim(tmp_para)[1]){
      files=tmp_para[i,.RES.COMBINE.COLS["PEP_SUMMARY_FILES"]]
      prjname=paste(tmp_para[i,c(.INPUT.COLS["SAMPLE"],.INPUT.COLS["BIOREP"])],
                    collapse="_",sep="")
      message("Processing file: ",prjname,"\n")
      logfile=paste(res$input_parameter$result.dir, "/", prjname, "-log.txt", 
                    collapse="", sep="");
      tmp_result=rbind(tmp_result,
          combineRun(pepFiles=files, fasta=res$input_parameter$fasta,
                     outPathFile=logfile, outdir=res$input_parameter$result.dir,
                     prefix=prjname))
    }
    
    tmp_result <- cbind(tmp_para,tmp_result)
    res$res_bioRep_level <- tmp_result
  }
  
  
  ## Combine result for the biological sample
  if(res$input_parameter$maxBioRep>=2){
    ## x is always data.frame
    tmp_para<-ddply(res$res_fraction_level,
                .variables=as.vector(c(.INPUT.COLS["SAMPLE"])),
                function(x){
                  y=data.frame(x=sum(x[,.RES.FRACTION.COLS["SPECTRUM_TOTAL"]]),
                               y=paste(x[,.RES.FRACTION.COLS["PEP_SUMMARY"]],
                                       sep="",collapse=";"))
                  colnames(y)=c(.RES.FRACTION.COLS["SPECTRUM_TOTAL"],
                                .RES.COMBINE.COLS["PEP_SUMMARY_FILES"])
                  y
                })
    tmp_result<-data.frame()
    for(i in 1:dim(tmp_para)[1]){
      files=tmp_para[i,.RES.COMBINE.COLS["PEP_SUMMARY_FILES"]]
      prjname=paste(tmp_para[i,c(.INPUT.COLS["SAMPLE"])],collapse="_",sep="")
      message("Processing file: ",prjname,"\n")
      logfile=paste(res$input_parameter$result.dir, "/", prjname, "-log.txt", 
                    collapse="", sep="");
      tmp_result=rbind(tmp_result,
          combineRun(pepFiles=files, fasta=res$input_parameter$fasta,
                     outPathFile=logfile, outdir=res$input_parameter$result.dir,
                     prefix=prjname))
    }
    
    tmp_result<-cbind(tmp_para,tmp_result)
    res$res_sample_level=tmp_result
  } else {
    ## if don't have bio replicate, the result for sample is the same with the 
    if(res$input_parameter$maxTechRep>=2){
      ## if don't have tech replicate
      res$res_sample_level = res$res_bioRep_level
    }else{
      if(res$input_parameter$maxFraction>=2){
        ## 
        res$res_sample_level = res$res_techRep_level
      }else{
        ## if don't have fracton
        
        res$res_sample_level = res$res_fraction_level
      }
    }
  }

  ## identification-independent metrics calculation
  ## if the spectra file is mz[X]ML format
  if(grepl(pattern="mz[x]*ml$",
           x=res$res_fraction_level[1,.INPUT.COLS["FILE"]],
           ignore.case=TRUE)){
    ms12.dir <- paste(res$input_parameter$result.dir,"/","ms12_metrics",sep="")
    dir.create(ms12.dir,recursive=TRUE,showWarnings=FALSE)
    res$ms12metrics <- calcMSQCMetrics(spectraList=spectralist,
                                     cpu=cpu,
                                     outdir=ms12.dir)
  
  }
  ## generate report
  ## search parameters: a table show the database search parameters
  ## TODO
  ## Identification summary table: a table show that the identification 
  ## result of all runs
  ## TODO
  ## make it an S3 for pretty printing
  class(res) <- c("list", "msQCres")
  saveRDS(res, file = file.path(outdir, "msQC.rds"))
  reportHTML(res)
  return(res)  
}

##' @title Load the result of \code{\link{msQCpipe}}
##' @description Load the result of \code{\link{msQCpipe}}
##' @param outdir The output directory of \code{\link{msQCpipe}}
##' @author Laurent Gatto \email{lg390@@cam.ac.uk}, 
##' Bo Wen \email{wenbo@@genomics.cn}
##' @examples 
##' zpqc <- system.file("extdata/qc.zip", package = "proteoQC")
##' unzip(zpqc)
##' qcres <- loadmsQCres("./qc")
loadmsQCres <- function(outdir) {
    x <- file.path(outdir, "msQC.rds")
    if (!file.exists(x))
        stop("Could not find ", x)
    res <- readRDS(x)
    return(res)
}


##' @title Print the information of msQCres object
##' @description Print the information of msQCres object
##' @param x A msQCres object
##' @param ... Additional parameters
##' @author Laurent Gatto \email{lg390@@cam.ac.uk}, 
##' Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' zpqc <- system.file("extdata/qc.zip", package = "proteoQC")
##' unzip(zpqc)
##' qcres <- loadmsQCres("./qc")
##' print.msQCres(qcres)
print.msQCres <- function(x, ...) {
    cat("An msQC results:\n")
    cat(" Results stored in", x[[1]]$outdir, "\n")
    cat(" Database:", x[[1]]$target.db, "\n")
    cat(" Run on", x[[1]]$date, "\n")
    n <- nrow(x[[2]])
    cat(" Design with", n, "samples:\n")
    if (n <= 3) {
        print(x[[2]][, 1:5])
    } else {
        print(x[[2]][1:2, 1:5])
        cat("  ... \n")
        print(x[[2]][n, 1:5])
    }
    cat("\nYou can now run reportHTML() to generate the QC report.\n")
}


##' @title Combine multiple results
##' @description Combine multiple results
##' @param pepFiles peptideSummary files
##' @param fasta database file
##' @param outPathFile out file
##' @param outdir output directory
##' @param prefix output prefix
##' @return A data.frame
##' @author Bo Wen \email{wenbo@@genomics.cn}
combineRun <- function(pepFiles,fasta,outPathFile,outdir,prefix){
    cmd_code <- paste("java -cp ",
                      paste('"',system.file("tool/tandemparser.jar", 
                                            package="proteoQC"),'"',
                            sep="",collapse=""), 
                      " cn.bgi.CombineRun ",
                      paste("\"",pepFiles,"\"",sep="",collapse=""), 
                      fasta,
                      outPathFile,
                      outdir,
                      prefix,
                      sep=" ", collapse=" ")
    ## cat(cmd_code,"\n")
    system(cmd_code)
    x<-read.table(outPathFile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    ## x<-as.data.frame(t(x))
    return(x)    
}


##' @title Run X!Tandem
##' @description Run X!Tandem
##' @title run X!Tandem
##' @description run X!Tandem
##' @param spectra MS/MS peak list file
##' @param fasta database file
##' @param outdir output directory
##' @param outprefix output file prefix
##' @param cpu The number of CPU used for X!Tandem
##' @param enzyme The ID of enzyme used for database searching. 
##' See \code{\link{showEnzyme}}.
##' @param xmx Set for parameter of "Java -Xmx".
##' @param varmod Variable modifications used for database searching.
##' See \code{\link{showMods}}.
##' @param fixmod Fixed modifications used for database searching.
##' See \code{\link{showMods}}.
##' @param tol The error window on experimental peptide mass values
##' @param tolu Units can be selected from: ppm, Daltons.
##' @param itol   Error window for MS/MS fragment ion mass values.
##' @param itolu  Units can be selected from: Daltons
##' @param miss Max miss clevage
##' @param refine Refine search, default is TRUE
##' @param ntt Default is 1
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @return a file path
runTandem=function(spectra="",fasta="",outdir="./",outprefix="",cpu=1,enzyme=1,
                    xmx=2,varmod=NULL,fixmod=NULL,refine=TRUE,ntt=1,
                   tol=10,tolu="ppm",itol=0.6,itolu="Daltons",miss=1){
  ##
  cat(format(Sys.time()),"\n")
  modsdb   = getMods()
  enzymedb = getEnzyme()
  cleavageSite = enzymedb$enzymestring[enzyme]
  varmods = ifelse(is.null(varmod),"",
                   paste(modsdb$modstring[varmod],collapse=","))
  fixmods = ifelse(is.null(fixmod),"",
                   paste(modsdb$modstring[fixmod],collapse=","))
  #taxonomy=rTTaxo(taxon="msqcfasta",format="peptide",URL=fasta)
  taxonomy=rTTaxo(taxon="msqcfasta",format="peptide",URL=fasta)
  outxmlname=paste(outdir,"/",basename(spectra),"_xtandem.xml",
                   collapse="",sep="")
  #taxonomy <- rTTaxo(
  #          taxon="yeast",
  #          format="peptide",
  #          URL=system.file("extdata/fasta/scd.fasta.pro", package="rTANDEM"))
  param <- rTParam()
  ##
  
  ##
  param <- setParamValue(param, 'list path', 'taxonomy information', taxonomy)
  param <- setParamValue(param, 'protein', 'taxon', value="msqcfasta")
  param <- setParamValue(param, 'list path', 'default parameters',
                         value=system.file("extdata/default_input.xml", 
                                           package="rTANDEM"))
  #param <- setParamValue(param, 'spectrum', 'path',
  #          value=system.file("extdata/test_spectra.mgf", package="rTANDEM"))
  param <- setParamValue(param, 'spectrum', 'path',value=spectra)
  param <- setParamValue(param, 'output', 'xsl path',
                         value=system.file("extdata/tandem-input-style.xsl", 
                                           package="rTANDEM"))
  param <- setParamValue(param, 'output', 'path',value=outxmlname)
  param <- setParamValue(param, 'output', 'maximum valid expectation value',
                         value=10)
  param <- setParamValue(param, 'output', 'parameters',value="yes")
  param <- setParamValue(param, 'output', 'results',value="all")
  param <- setParamValue(param, 'output', 'path hashing',value="no")
  
  param <- setParamValue(param,"spectrum","fragment monoisotopic mass error",
                         value=itol)
  ##The value for this parameter may be 'Daltons' or 'ppm': all other values are
  ##ignored
  param<-setParamValue(param,"spectrum",
                       "fragment monoisotopic mass error units",
                         value=itolu)
  param <-setParamValue(param,"spectrum","parent monoisotopic mass error plus",
                         value=tol)
  param <-setParamValue(param,"spectrum","parent monoisotopic mass error minus",
                         value=tol)
  param <-setParamValue(param,"spectrum","parent monoisotopic mass error units",
                         value=tolu)
  ##The value for this parameter may be 'Daltons' or 'ppm': all other values 
  ##are ignored
  param <- setParamValue(param,"spectrum","maximum parent charge",value=8)
  ##don't aotu-decoy search, just search against a database which contained 
  ##decoy sequences
  param <- setParamValue(param,"scoring","include reverse",value="no")
  param <- setParamValue(param,"scoring","maximum missed cleavage sites",
                         value=miss)
  param <- setParamValue(param,"spectrum","threads",value=cpu)
  
  param <- setParamValue(param,"protein","cleavage site",value=cleavageSite)
  param <- setParamValue(param,"residue","potential modification mass",
                         value=varmods)
  param <- setParamValue(param,"residue","modification mass",value=fixmods)
  
  param <- setParamValue(param,"refine",value=ifelse(refine,"yes","no"))
  
  if(refine){
    param <- setParamValue(param,"refine","unanticipated cleavage",value=ifelse(ntt==2,"no","yes"))
  }
   
  
  
  result.path <- tandem(param)
  
  ## use the java parser to extract the information
  ##the information output file: log.txt
  logfile = paste(outdir,"/",outprefix,"-log.txt",collapse="",sep="");
  
  tandemparser=paste("java ",paste("-Xmx",xmx,"G",sep="")," -jar ",
                     paste('"',system.file("tool/tandemparser.jar", 
                                           package="proteoQC"),'"',
                           sep="",collapse=""), 
                     result.path, fasta ,outprefix,outdir,logfile,
                     collapse=" ",sep=" ")
  outfile=system(command=tandemparser,intern=TRUE)
  cat(outfile,"\n")
  x<-read.table(logfile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  colnames(x)<-c(.RES.FRACTION.COLS["SPECTRUM_TOTAL"], 
                 .RES.FRACTION.COLS["SPECTRUM"],
                 .RES.FRACTION.COLS["PEPTIDE"],
                 .RES.FRACTION.COLS["PROTEIN"],
                 #.RES.FRACTION.COLS["MSMS_ERROR"],
                 .RES.FRACTION.COLS["PEP_SUMMARY"],
                 .RES.FRACTION.COLS["PRO_SUMMARY"])
  
  cat(format(Sys.time()),"\n")
  return(x)
}


##' @title Shown all modifications 
##' @description Shown all modifications 
##' @return A data frame which contains all of the modifications
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' showMods()
showMods=function(){
  modsdb = read.delim(system.file("config/mods.txt", package="proteoQC"),
      header=TRUE,sep="\t",stringsAsFactors=FALSE)
  modsdb
  
}

##' @title Shown all enzymes
##' @description Shown all enzymes
##' @return A data frame which contains all of the enzymes
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' showEnzyme()
showEnzyme=function(){
  enzymedb = read.delim(system.file("config/enzyme.txt",
      package="proteoQC"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  enzymedb
  
}

##' @title Get the modification list
##' @description Get the modification list
##' @return A data frame which contains all of the modifications
##' @author Bo Wen \email{wenbo@@genomics.cn}
getMods=function(){
  modsdb = read.delim(system.file("config/mods.txt", package="proteoQC"),
      header=TRUE,sep="\t",stringsAsFactors=FALSE)
  return(modsdb)
}

##' @title Get the enzymes list
##' @description Get the enzymes list
##' @return A data frame which contains all of the enzymes
##' @author Bo Wen \email{wenbo@@genomics.cn}
getEnzyme=function(){
  enzymedb = read.delim(system.file("config/enzyme.txt",
      package="proteoQC"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  return(enzymedb)
}


##' @title Create target-decoy database
##' @description Create target-decoy database
##' @param fa target database
##' @param outdb output target-decoy database
##' @return target-decoy database file name
##' @author Bo Wen \email{wenbo@@genomics.cn}
createTargetDecoyDB=function(fa,outdb){
  
  seq.db  <- read.fasta(file=fa,seqtype="AA")
 
  seq.cnt <- read.fasta(file=system.file("db/crap.fasta",package="proteoQC"),
                        seqtype="AA")
  ## The id and description may be separated by "\t" 
  seq.db <- lapply(seq.db,function(x){
    id <- attr(x,which="name")
    attr(x,which="name")<-gsub(pattern="\\s+.*$",replacement="",x=id,perl=TRUE)
    x
  })
  seq.cnt <- lapply(seq.cnt,function(x){
    id <- attr(x,which="name")
    attr(x,which="name")<-gsub(pattern="\\s+.*$",replacement="",x=id,perl=TRUE)
    attr(x,which="name")<-paste("cnt_msQC_",attr(x,which="name"),sep="")
    x
  })
  
  ## Get the reverse seq.
  seq.revdb <- lapply(seq.db,function(x){
    id <- attr(x,which="name")
    attr(x,which="name") <- paste0("###REV###",id)
    y <- rev(x)
    attributes(y) <- attributes(x)
    y
  })
  seq.revcnt <- lapply(seq.cnt,function(x){
    id <- attr(x,which="name")
    attr(x,which="name") <- paste0("###REV###",id)
    y <- rev(x)
    attributes(y) <- attributes(x)
    y
  })
  
  write.fasta(sequences=seq.db,names=sapply(seq.db,attr,which="name"),
              file.out=outdb)
  write.fasta(sequences=seq.cnt,names=sapply(seq.cnt,attr,which="name"),
              file.out=outdb,open="a")
  write.fasta(sequences=seq.revdb,names=sapply(seq.revdb,attr,which="name"),
              file.out=outdb,open="a")
  write.fasta(sequences=seq.revcnt,names=sapply(seq.revcnt,attr,which="name"),
              file.out=outdb,open="a")
  
  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Definitions.
###

# Definition of column names corresponding to input experiment files
.INPUT.COLS <- c(FILE="file",
                 SAMPLE="sample",
                 BIOREP="bioRep",
                 TECHREP="techRep",
                 FRACTION="fraction")

.RES.FRACTION.COLS <- c(SPECTRUM_TOTAL="spectrum_total",
                        SPECTRUM="spectrum",
                        PEPTIDE="peptide",
                        PROTEIN="protein",
                        MSMS_ERROR="msms_error",
                        PEP_SUMMARY="peptide_summary",
                        PRO_SUMMARY="protein_summary")
.RES.COMBINE.COLS <- c(PEP_SUMMARY_FILES="pep_summary_files")

.PepSummaryColClass <- c("character", #index
                         "numeric", #evalue
                         "numeric", #charge
                         "numeric", #mass
                         "numeric", #mz
                         "numeric", #delta_da
                         "numeric", #delta_ppm
                         "character", #peptide
                         "numeric", #isdecoy
                         "numeric", #miss
                         "character", #protein
                         "character", #rt
                         "character", #mods
                         "numeric", #ms2peaks
                         #"numeric", #ms2delta
                         "character", #ms2delta
                         "numeric", #FDR
                         "numeric"  #Qvalue
                         )
.ProSummaryColClass <- c("character", #Accession
                         "numeric", #Mass
                         "character", #PeptideSeqs
                         "character", #PepIsUniqe
                         "numeric", #NumOfUniqPeps
                         "numeric", #NumOfUniqSpectra
                         "character", #Spectra
                         "character", #SameSet
                         "character" #Description
                        )



