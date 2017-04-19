library("msQC")

trimMgf <- function(f, m = 1/10, overwrite = FALSE) {
    message("Reading ", f)
    x <- readLines(f)
    beg <- grep("BEGIN IONS", x)
    end <- grep("END IONS", x)
    n <- length(beg)
    message("Sub-setting to ", m)
    i <- sort(sample(n, floor(n * m)))
    k <- unlist(mapply(seq, from = beg[i], to = end[i]))
    if (overwrite) {
        unlink(f)
        message("Writing ", f)
        writeLines(x[k], con = f)
        return(f)
    } else {
        g <- sub(".mgf", "_small.mgf", f)
        message("Writing ", g)
        writeLines(x[k], con = g)
        return(g)
    }    
}

set.seed(1)

trimMgf("Thermo_Hela_PRTC_1.mgf")
trimMgf("Thermo_Hela_PRTC_2.mgf")
trimMgf("Thermo_Hela_PRTC_3.mgf")

qcres <- msQCpipe(spectralist = "design_small.txt",
                  fasta = "swissprot_human_canonical_19_09_12.fasta",
                  outdir = "./qc",
                  miss  = 0,
                  enzyme = 1, varmod = 2, fixmod = 1,
                  tol = 10, itol = 0.6, cpu = 2,
                  mode = "identification")

unlink("./qc/database/target_decoy.fasta")
unlink("./qc/result/*_xtandem.xml")

## reportHTML("./qc")
## unlink("./qc/report/*")
## unlink("./qc/qc_report.*")

zip("qc.zip", "qc")
