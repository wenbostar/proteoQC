
test_loadmsQCres=function(){
  zpqc <- system.file("extdata/qc.zip", package = "proteoQC")
  unzip(zpqc)
  qcres <- loadmsQCres("./qc")
  checkEquals(dim(qcres$res_fraction_level)[1],8)
}