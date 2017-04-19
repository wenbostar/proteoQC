[![Bioconductor devel build Status](http://bioconductor.org/shields/build/devel/bioc/proteoQC.svg)](http://bioconductor.org/packages/devel/bioc/html/proteoQC.html) 

```{r env, echo=FALSE,}
suppressPackageStartupMessages(library("BiocStyle"))
knitr::opts_chunk$set(comment=NA)
```

# `proteoQC`
proteoQC is an R package for proteomics data quality assessment. This package creates an HTML format QC report for MS/MS-based proteomics data. The report is intended to allow the user to quickly assess the quality of proteomics data.

## Installation

To install `r Biocpkg("proteoQC")`

```{r install, eval = FALSE}
library("BiocInstaller")
biocLite("proteoQC")
```

If you need the github version (not recommended unless you know what
you are doing)

```{r installgh, eval = FALSE}
biocLite("wenbostar/proteoQC")
```


## Contribution

Contributions to the package are more than welcome. 
