setMethod("genotype", "SnpQSet",
          function(object, correction=NULL, recalibrate=TRUE,
                   minLLRforCalls=c(5, 1, 5), verbose=TRUE,
                   correctionFile=NULL, prefix="tmp.crlmm.",
                   balance=1.5)

          crlmm(object, correction=correction,
                recalibrate=recalibrate,
                minLLRforCalls=minLLRforCalls, verbose=verbose,
                correctionFile=correctionFile, prefix=prefix,
                balance=balance)
        )

setMethod("genotype", "character",
          function(object, tmpdir, batch_size=40000,
                      balance=1.5, minLLRforCalls=c(5, 1, 5),
                      recalibrate=TRUE, verbose=TRUE, pkgname)

          genotypeOne(object, tmpdir=tmpdir, batch_size=batch_size,
                      balance=balance, minLLRforCalls=minLLRforCalls,
                      recalibrate=recalibrate, verbose=verbose,
                      pkgname=pkgname)
          
        )
