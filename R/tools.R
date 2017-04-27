
##' @title Protein inference
##' @description Protein inference
##' @param file A file containing the information of peptides to proteins.
##' @param db A protein database of fasta format.
##' @param pepColName The column name of peptide sequence.
##' @param proColName The column name of protein ID.
##' @param spectrumColName The column name of spectrum index.
##' @param proSep The separator of protein ID, default is "".
##' @param outfile The output file name of protein group result.
##' @param xmx JAVA -Xm
##' @return NULL
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' pep.zip <- system.file("extdata/pep.zip", package = "proteoQC")
##' unzip(pep.zip)
##' proteinGroup(file = "pep.txt", outfile = "pg.txt")
proteinGroup=function(file=NULL,db="",pepColName="peptide",
                      proColName="protein",spectrumColName="index",
                      proSep=";",outfile=NULL,xmx=1){
  ## get the java tool
  pgbin <- paste("java ",paste("-Xmx",xmx,"G",sep="")," -cp ",
                     paste('"',system.file("tool/tandemparser.jar", 
                                           package="proteoQC"),'"',
                           sep="",collapse=""), " cn.bgi.ProteinGroup ",
                     file, outfile ,paste('"',db,'"',sep="",collapse = ""),
                     pepColName, spectrumColName, proColName,
                     paste('"',proSep,'"',sep="",collapse = ""),
                     collapse=" ",sep=" ")
  ## run
  outfile=system(command=pgbin,intern=TRUE)
  cat(outfile,"\n")    
}



##' @title Charge distribution
##' @description Read the charge information from mgf file
##' @param mgf A file of mgf.
##' @return A data.frame object
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' mgf.zip <- system.file("extdata/mgf.zip", package = "proteoQC")
##' unzip(mgf.zip)
##' charge <- chargeStat("test.mgf")
chargeStat=function(mgf=NULL){
  msmsdata <- readMgfData(mgf)
  charge_result <- precursorCharge(msmsdata)
  dat <- as.data.frame(table(charge_result))
  names(dat) <- c("Charge","Number")
  return(dat);
}

##' @title Calculate the labeling efficiency of isobaric labeling data 
##' @description Calculate the labeling efficiency of isobaric labeling data
##' @param ms MS/MS file.
##' @param reporter Isobaric tag class, 1=iTRAQ-4plex, 2=iTRAQ-8plex, 3=TMT-6plex.
##' 4=TMT-10plex.
##' @param plot Logical value
##' @return A vector object
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' mgf.zip <- system.file("extdata/mgf.zip", package = "proteoQC")
##' unzip(mgf.zip)
##' a <- labelRatio("test.mgf",reporter=2)
labelRatio=function(ms=NULL,reporter=1,plot=TRUE){
  msmsdata <- readMgfData(ms)
  if(reporter==1){
      reporterClass <- iTRAQ4
  }else if(reporter==2){
      reporterClass <- iTRAQ8
  }else if(reporter==3){
      reporterClass <- TMT6
  }else if(reporter==4){
      reporterClass <- TMT10
  }
  result <- as.data.frame(exprs(quantify(msmsdata, reporters = reporterClass, method = "max")))

  if(plot){
      
    dat <- gather(result,"Tag","Intensity") 
    nspectra <- nrow(result)
    databar <- group_by(dat,Tag) %>% dplyr::summarise(n=sum(!is.na(Intensity)))
    databar$ratio <- databar$n/nspectra
    databar$label <- sprintf("%.2f%%",databar$ratio*100)
    
    p <- ggplot(data = databar, aes(x=Tag,y=ratio, fill=Tag))+
      geom_bar(stat="identity",width=.5)+
      xlab("Isobaric Tag")+
      ylab("Ratio of labeling")+
      theme(legend.position = "none")+
      geom_text(mapping = aes(label=label),vjust=-0.8,size=3)
    print(p)
    
    p <- ggplot(data = dat, aes(x=Tag,y=log2(Intensity), fill=Tag))+
        geom_boxplot(width=0.4)+
        xlab("Isobaric Tag")+
        ylab("log2(Intensity)")+
        theme(legend.position = "none")
        
    print(p)
  }
  
  return(result);
}
