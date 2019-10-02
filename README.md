[![Build Status](https://travis-ci.org/wenbostar/proteoQC.svg?branch=master)](https://travis-ci.org/wenbostar/proteoQC) 
[![Bioconductor release build Status](http://bioconductor.org/shields/build/release/bioc/proteoQC.svg)](http://bioconductor.org/packages/release/bioc/html/proteoQC.html) 
[![Bioconductor devel build Status](http://bioconductor.org/shields/build/devel/bioc/proteoQC.svg)](http://bioconductor.org/packages/devel/bioc/html/proteoQC.html) 


# `proteoQC`
*[proteoQC](http://bioconductor.org/packages/proteoQC)* is an R package for proteomics data quality assessment. This package creates an HTML format QC report for MS/MS-based proteomics data. The report is intended to allow the user to quickly assess the quality of proteomics data. The official page is the Bioconductor landing page
([release](http://www.bioconductor.org/packages/release/bioc/html/proteoQC.html)
or
[devel](http://www.bioconductor.org/packages/devel/bioc/html/proteoQC.html)
versions).

## Installation

To install *[proteoQC](http://bioconductor.org/packages/proteoQC)*


```{r install, eval = FALSE}
library("BiocInstaller")
biocLite("proteoQC")
```

If you need the github version (not recommended unless you know what
you are doing)

```{r installgh, eval = FALSE}
biocLite("wenbostar/proteoQC")
```
## Citation

To cite the `proteoQC` package in publications, please use:

> Wen B and Gatto L (2017). proteoQC: An R package for proteomics data quality control. R package version 1.15.0, https://github.com/wenbostar/proteoQC.

## List of citations

`proteoQC` has been cited in the following manuscripts:
1. Gatto, Laurent, et al. "Visualization of proteomics data using R and Bioconductor." Proteomics 15.8 (2015): 1375-1389.
2. Bittremieux, Wout, et al. "Computational quality control tools for mass spectrometry proteomics." Proteomics 17.3-4 (2017).
3. Samandi, Sondos, et al. "Deep transcriptome annotation enables the discovery and functional characterization of cryptic small proteins." eLife 6 (2017): e27860.
4. Belghit I, Lock E J, Fumière O, et al. Species-Specific Discrimination of Insect Meals for Aquafeeds by Direct Comparison of Tandem Mass Spectra. Animals, 2019, 9(5): 222.
5. Walzer M., Vizcaíno J.A. (2020) Review of Issues and Solutions to Data Analysis Reproducibility and Data Quality in Clinical Proteomics. In: Matthiesen R. (eds) Mass Spectrometry Data Analysis in Proteomics. Methods in Molecular Biology, vol 2051. Humana, New York, NY


## Contribution

Contributions to the package are more than welcome. 
