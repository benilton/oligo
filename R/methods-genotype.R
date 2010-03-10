snpArrays1 <- c("pd.mapping50k.xba240", "pd.mapping50k.hind240",
                "pd.mapping250k.nsp", "pd.mapping250k.sty")
snpArrays2 <- c("pd.genomewidesnp.5", "pd.genomewidesnp.6")

crlmm <- function(filenames, outdir, batch_size=40000, balance=1.5,
                    minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                    verbose=TRUE, pkgname, reference=TRUE){
  stopifnot(!(missing(filenames) | missing(outdir)))
  if(missing(pkgname)) pkgname <- cleanPlatformName(readCelHeader(filenames[1])$chiptype)
  requireAnnotation(pkgname, verbose=verbose)

  
  if(pkgname %in% snpArrays1){
    justCRLMMv3(filenames, outdir, batch_size=batch_size,
              recalibrate=recalibrate, minLLRforCalls=minLLRforCalls,
              balance=balance, verbose=verbose)
  }else if (pkgname %in% snpArrays2){
    genotypeOne(filenames, outdir, batch_size=batch_size,
                balance=balance, minLLRforCalls=minLLRforCalls,
                recalibrate=recalibrate, verbose=verbose,
                pkgname=pkgname, reference=reference)
  }else{
    stop("CRLMM not implemented for ", pkgname)
  }
}
