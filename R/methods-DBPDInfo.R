setMethod("getX", "DBPDInfo",
          function(object, type){
            stopifnot(!missing(type))
            sql <- "SELECT fid, x FROM"
            if (type == "pm"){
              sql <- paste(sql, "pmfeature")
            }else if(type == "bg"){
              sql <- paste(sql, "bgfeature")
            }else{
              stop("Method not implemented for type ", type)
            }
            res <- dbGetQuery(db(object), sql)
            res[order(res[["fid"]]), "x"]
          })

setMethod("getY", "DBPDInfo",
          function(object, type){
            stopifnot(!missing(type))
            sql <- "SELECT fid, y FROM"
            if (type == "pm"){
              sql <- paste(sql, "pmfeature")
            }else if(type == "bg"){
              sql <- paste(sql, "bgfeature")
            }else{
              stop("Method not implemented for type ", type)
            }
            res <- dbGetQuery(db(object), sql)
            res[order(res[["fid"]]), "y"]
          })

## from oligoClasses
setMethod("pmindex", "DBPDInfo",
          function(object, subset=NULL) { 
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            tmp <- dbGetQuery(db(object),
                              "select fid from pmfeature")[[1]]
            sort(tmp)
          })

setMethod("mmindex", "DBPDInfo",
          function(object, subset=NULL){
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            tmp <- dbGetQuery(db(object),
                       "select fid from mmfeature")[[1]]
            sort(tmp)
          })

setMethod("probeNames", "DBPDInfo",
          function(object, subset=NULL) {
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            sql <- "select man_fsetid, fid from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            tmp[order(tmp[["fid"]], tmp[["man_fsetid"]]), "man_fsetid"]
          })

setMethod("geneNames", "DBPDInfo",
          function(object) {
              unique(probeNames(object))
          })

setMethod("pmSequence", "AffySNPPDInfo",
          function(object){
            sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("mmSequence", "AffySNPPDInfo",
          function(object){
            sql <- "select seq from sequence, mmfeature where mmfeature.fid=sequence.fid order by mmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("bgSequence", "DBPDInfo",
          function(object){
            theFile <- file.path(system.file(package=annotation(object)), "data", "bgSequence.rda")
            if (file.exists(theFile)){
              load(theFile)
              return(bgSequence[["sequence"]])
            }else{
              warning("bgSequence.rda file is not available for this pdInfo pkg.")
              return(NULL)
            }
          })

setMethod("pmSequence", "DBPDInfo",
          function(object){
            theFile <- file.path(system.file(package=annotation(object)), "data", "pmSequence.rda")
            if (file.exists(theFile)){
              load(theFile)
              return(pmSequence[["sequence"]])
            }else{
              warning("pmSequence.rda file is not available for this pdInfo pkg.")
              return(NULL)
            }
          })

setMethod("mmSequence", "DBPDInfo",
          function(object){
            theFile <- file.path(system.file(package=annotation(object)), "data", "mmSequence.rda")
            if (file.exists(theFile)){
              load(theFile)
              return(mmSequence[["sequence"]])
            }else{
              warning("mmSequence.rda file is not available for this pdInfo pkg.")
              return(NULL)
            }
          })

setMethod("pmOffset", "AffySNPPDInfo",
          function(object){
            sql <- "select offset from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmFragmentLength", "AffySNPPDInfo",
          function(object){
            sql <- "select fid, fragment_length from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            idx <- order(tmp[["fid"]])
            tmp[idx, "fragment_length"]
          })

setMethod("pmAllele", "AffySNPPDInfo",
          function(object){
            sql <- "select allele from pmfeature order by fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmStrand", "AffySNPPDInfo",
          function(object){
            sql <- "select strand from pmfeature order by fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmPosition", "ExpressionPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmPosition", "TilingPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("bgindex", "DBPDInfo",
          function(object, subset=NULL){
            if (!is.null(subset))
              warning("Subset not implemented (yet). Returning everything.")
            sql <- "SELECT fid FROM bgfeature"
            tmp <- dbGetQuery(db(object), sql)[[1]]
            sort(tmp)
          })
