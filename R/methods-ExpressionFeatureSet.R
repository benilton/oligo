setMethod("backgroundCorrect", "ExpressionFeatureSet",
          function(object, method="rma", copy=TRUE, verbose=TRUE, ...){
              method <- match.arg(method, c("rma", "mas"))
              if (verbose) message("Background correcting... ",
                                   appendLF=FALSE)
              if (method == "rma"){
                  out <- callNextMethod(object=object, method="rma",
                                        copy=TRUE, verbose=FALSE)
              }else if (method == "mas"){
                  dots <- list(...)
                  if (is.null(dots[["griddim"]])){
                      griddim <- 16
                  }else{
                      griddim <- dots[["griddim"]]
                  }
                  out <- bgMAS(object, griddim=griddim)
              }
              if (verbose) message("OK")
              return(out)
          })

bgMAS <- function(object, griddim=16){
    ## mod by BC
   nchips <- ncol(object)

   pm.index <- pmindex(object)
   mm.index <- mmindex(object)

   ## some chips have some probesets without MM probes
   ## which will return an NA in mm.index

   mm.index <- mm.index[!is.na(mm.index)]

   ## rows/cols of the *chip*
   rows <- geometry(object)[1]
   cols <- geometry(object)[2]

   allintensities <- exprs(object)[c(pm.index, mm.index), ]

   # note that the indexing is +1 more than you'd expect because
   # the c code expects it that way
   ## (note about the remark above: R indexing starts at 1 and not at 0,
   ## that's why the indexing is done this way. The package is primarily done to
   ## be used with R...)

   allx <- c(pm.index-1, mm.index-1) %% rows +1
   ally <- c(pm.index-1, mm.index-1) %/% rows + 1

   nprobes <- length(allx)

   ## affy_background_adjust_R could be in preprocessCore
   corrected <- matrix(.C("affy_background_adjust_R",
                          as.double(as.vector(allintensities)), as.integer(allx), as.integer(ally),
                          as.integer(nprobes), as.integer(nchips), as.integer(rows), as.integer(cols),
                          as.integer(griddim), PACKAGE="oligo")[[1]],
                       nprobes, nchips)

   ## modified by BC
   newExprs <- exprs(object)
   newExprs[c(pm.index, mm.index), ] <- corrected
   out <- object
   exprs(out) <- newExprs

   ## and what with the 'non pm or mm' probes ?
   ## answer: they are not used per Affymetrix Statistical Algorithms Description Document.

   return(out)

 }

setMethod("rma", "ExpressionFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            ## get pmi and pnVec
            pmi <- pmindex(object)
            pnVec <- probeNames(object)
            tbls <- dbListTables(db(object))
            if (manufacturer(object) == "Affymetrix" && "bgfeature" %in% tbls){
              sql <- paste("SELECT man_fsetid, fid",
                           "FROM bgfeature",
                           "INNER JOIN featureSet",
                           "USING(fsetid)")
              tmpQcPm <- dbGetQuery(db(object), sql)
              pmi <- c(pmi, tmpQcPm[["fid"]])
              pnVec <- c(pnVec, tmpQcPm[["man_fsetid"]])
            }
            idx <- order(pnVec)
            pnVec <- pnVec[idx]
            pmi <- pmi[idx]
            rm(idx)
            
            theClass <- class(exprs(object))

            if ("matrix" %in% theClass){
              pms <- exprs(object[pmi,])
              dimnames(pms) <- NULL
              colnames(pms) <- sampleNames(object)
              exprs <- basicRMA(pms, pnVec, normalize, background)
            }else if ("ff_matrix" %in% theClass){
              pms <- ffSubset(rows=pmi, object=exprs(object), prefix="pm-")
              exprs <- basicRMAbo(pms, pnVec, background=background, normalize=normalize)
              finalizer(pms) <- "delete"
              rm(pms)
            }else{
              stop("basicRMA not implemented for '", theClass, "' objects.")
            }
            
            out <- new("ExpressionSet")
            slot(out, "assayData") <- assayDataNew(exprs=exprs)
            slot(out, "phenoData") <- phenoData(object)
            slot(out, "featureData") <- basicAnnotatedDataFrame(exprs, byrow=TRUE)
            slot(out, "protocolData") <- protocolData(object)
            slot(out, "annotation") <- slot(object, "annotation")
            if (validObject(out)){
              return(out)
            }else{
              stop("Resulting object is invalid.")
            }
          })


setMethod("getNetAffx", "ExpressionSet",
          function(object, type="probeset"){
              type <- match.arg(type, c("probeset", "transcript"))
              fname <- ifelse(type == "probeset",
                              "Probeset",
                              "Transcript")
              fname <- paste("netaffx", fname, ".rda", sep="")
              fname <- file.path(system.file("extdata", package=annotation(object)),
                                 fname)

              if (!file.exists(fname))
                  stop("NetAffx Annotation not available in '",
                       annotation(object), "'. Consider using 'biomaRt'.")

              obj <- load(fname)
              fdata <- get(obj)
              rm(list=obj)
              fns <- featureNames(object)
              theSet <- fdata[fns,]
              featureNames(theSet) <- fns
              theSet
          })

