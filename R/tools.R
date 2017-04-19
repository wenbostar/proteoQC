
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
##' @export
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
##' @return A vector object
##' @export
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' mgf.zip <- system.file("extdata/mgf.zip", package = "proteoQC")
##' unzip(mgf.zip)
##' charge <- chargeStat("test.mgf")
chargeStat=function(mgf=NULL){
  result <- .Call('ChargeCount_Cpp', PACKAGE = 'proteoQC', mgf)
  charge <- unlist(result)
  names(charge) <- gsub(pattern = "\\+.*$",replacement = "",x=names(charge))
  return(charge);
}

##' @title Calculate the labeling efficiency of isobaric labeling data 
##' @description Calculate the labeling efficiency of isobaric labeling data
##' @param ms MS/MS file.
##' @param iClass Isobaric tag class, 1=iTRAQ-8plex.
##' @param delta The mass error for reporter matching.
##' @param plot Logical value
##' @return A vector object
##' @export
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' mgf.zip <- system.file("extdata/mgf.zip", package = "proteoQC")
##' unzip(mgf.zip)
##' a <- labelRatio("test.mgf")
labelRatio=function(ms=NULL,iClass=1,delta=0.05,plot=TRUE){
  result <- .Call('LableRatio_Cpp', PACKAGE = 'proteoQC', ms,iClass,delta)
  result <- unlist(result)
  
  if(plot){
    label.names <- NULL
    if(iClass==1){
      label.names <- paste("I",c(113:119,121),sep="")
    }
    label.names <- c(label.names,"none")
    dat <- sapply(label.names,function(x){
      sum(result[grep(pattern = x,x = names(result))])})
    dat <- data.frame(x=dat/sum(result),label=names(dat))
    dat$ratio <- sprintf("%.2f%%",dat$x*100)
    p <- ggplot(data = dat, aes(x=label,y=x, fill=label))+
      geom_bar(stat="identity",width=.5)+
      xlab("Isobaric Tag")+
      ylab("Ratio of labeling")+
      theme(legend.position = "none")+
      geom_text(mapping = aes(label=ratio),vjust=-0.8,size=5)
    print(p)
  }
  
  return(result);
}
