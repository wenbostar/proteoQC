library("msQC")
library("R.utils")

library("rpx")
px <- PXDataset("PXD000864")

mgfs <- grep("mgf", pxfiles(px), value = TRUE)
mgfs16 <- grep("-0[1-6]-", mgfs, value=TRUE)

mgffiles <- pxget(px, mgfs16)
mgffiles <- sapply(mgffiles, gunzip)

fas <- pxget(px, "TTE2010.zip")
fas <- unzip(fas)
fas

sample <- 1:length(mgfs16)
techrep <- rep(1:3, 12)
biorep <- rep(1, 36)
frac <- rep((rep(1:6, each = 3)), 2)
des <- cbind(mgf = sub(".gz", "", mgfs16),
             sample = sample,
             bioRep = biorep, techRep = techrep,
             fraction = frac)
write.table(des, sep = " ", row.names=FALSE,
            quote = FALSE,
            file = "./PXD000864-design.txt")

qcres <- msQCpipe(spectralist = design,
                  fasta = fas, 
                  outdir = "./qc",
                  miss  = 0,
                  enzyme = 1, varmod = 2, fixmod = 1,
                  tol = 10, itol = 0.6, cpu = 2,
                  mode = "identification")

zip("./qc_36.zip", "./qc")

