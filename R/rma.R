######################################################
#
# rma - RMA interface to c code
#
# the RMA method implemented in c code
#
# this code serves as interface to the c code.
# currently
# implemented (version 0.25) background correction
#
# Background correction code has been added.
#
# note this function does not leave the supplied
# AffyBatch unchanged if you select DESTRUCTIVE=TRUE. this is
# for memory purposes but can be quite
# dangerous if you are not careful. Use destructive=FALSE if this is
# deemed likely to be a problem.
#
# UPDATE: note that the affybatch is now not affected if you use
# destructive=TRUE and you might actually save a little memory.
# the destructive refers only to Plobs, which would be destroyed.
#
# History
#
# Feb 22, 2004 - activated subset. In is now possible to
#                do the entire RMA procedure using a subset of probesets
#
#
# May 20, 2005 - adapted for oligo package [by Rafael Irizarry]
#
#
#
########################################################

##subset does now work yet
rma <- function(object, subset=NULL, verbose=TRUE, destructive = TRUE,
                normalize=TRUE, background=TRUE, bgversion=2, ...)
{
    pnms <- probeNames(object, subset)
    rows <- length(pnms)
    cols <- length(object)

    if (is.null(subset)){
        ngenes <- length(geneNames(object))
    } else {
        ngenes <- length(subset)
    }

    ## background correction
    bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

    fsetPMs <- pm(object, subset)
    if (destructive){
        exprs <- .Call("rma_c_complete", fsetPMs, fsetPMs,
                       pnms, ngenes, body(bg.dens),
                       new.env(), normalize, background, bgversion, PACKAGE="oligo")
    } else {
        exprs <- .Call("rma_c_complete_copy", fsetPMs,
                       fsetPMs, pnms, ngenes,
                       body(bg.dens), new.env(), normalize, background, bgversion,
                       PACKAGE="oligo")
    }
    colnames(exprs) <- sampleNames(object)
    ## to be fixed later, besides which don't believe much in nominal
    ## se's with medianpolish
    se.exprs <- array(NA, dim(exprs))
    dimnames(se.exprs) <- dimnames(exprs)
    ## NOTE: when creating a new ExpressionSet in parts as below, _order
    ## matters_.  Each replacement function induces a copy.  Hence, it
    ## is best to put the larger data items in last.
    out <- new("ExpressionSet")
    annotation(out) <- annotation(object)
    experimentData(out) <- experimentData(object)
    phenoData(out) <- phenoData(object)
    assayData(out) <- assayDataNew(exprs=exprs, se.exprs=se.exprs)
    ## FIXME: should we call validObject here?
    return(out)
}
