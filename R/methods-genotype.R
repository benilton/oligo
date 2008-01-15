crlmm <- function(filenames, outdir, batch_size=40000, balance=1.5,
                    minLLRforCalls=c(5, 1, 5), recalibrate=TRUE,
                    verbose=TRUE, pkgname){
  stopifnot(!(missing(filenames) | missing(outdir)))
  if(missing(pkgname)) pkgname <- cleanPlatformName(readCelHeader(filenames[1])$chiptype)
  require(pkgname, character.only=TRUE)
  if(pkgname %in% c("pd.mapping50k.xba240", "pd.mapping50k.hind240",
                    "pd.mapping250k.nsp", "pd.mapping250k.sty")){
    justCRLMMv2(filenames, outdir, batch_size=batch_size,
              recalibrate=recalibrate, minLLRforCalls=minLLRforCalls,
              balance=balance, verbose=verbose)
  }else if (pkgname %in% c("pd.genomewidesnp.5", "pd.genomewidesnp.6")){
    genotypeOne(filenames, outdir, batch_size=batch_size,
                balance=balance, minLLRforCalls=minLLRforCalls,
                recalibrate=recalibrate, verbose=verbose,
                pkgname=pkgname)
  }
}
