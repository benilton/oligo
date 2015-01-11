getFSetInfo <- function(obj, id, fields){
    conn <- db(obj)
    tbl <- paste0('featureSet', id)
    if (! (tbl %in% dbListTables(conn))){
        avail <- paste(dbListTables(conn), collapse=', ')
        message('Available tables: ', avail)
        stop('Table ', tbl, ' does not exist.')
    }
    if (!is.null(fields)){
        fields <- unique(c('fsetid', 'man_fsetid', fields))
        fields <- paste0(fields, collapse=', ')
    }else{
        fields <- '*'
    }
    sql <- paste('SELECT', fields, 'FROM', tbl)
    dbGetQuery(conn, sql)
}

getMPSInfo <- function(obj, id, fields, type='pm'){
    conn <- db(obj)
    tbl <- paste0('mps', id, type)
    fsettbl <- paste0('featureSet', id)
    if (! (tbl %in% dbListTables(conn))){
        avail <- paste(dbListTables(conn), collapse=', ')
        message('Available tables: ', avail)
        stop('Table ', tbl, ' does not exist.')
    }
    if (!missing(fields)){
        fields <- unique(c('fsetid', 'man_fsetid', 'fid', fields))
        fields <- paste0(fields, collapse=', ')
    }else{
        fields <- '*'
    }
    sql <- paste('SELECT', fields, 'FROM', tbl, 'INNER JOIN', fsettbl, 'USING(fsetid)')
    dbGetQuery(conn, sql)
}

getFeatureInfo <- function(obj, fields, type='pm'){
    conn <- db(obj)
    tbl <- paste0(type, 'feature')
    if (! (tbl %in% dbListTables(conn))){
        avail <- paste(dbListTables(conn), collapse=', ')
        message('Available tables: ', avail)
        stop('Table ', tbl, ' does not exist.')
    }
    if (!missing(fields)){
        fields <- unique(c('fid', fields))
        fields <- paste0(fields, collapse=', ')
    }else{
        fields <- '*'
    }
    sql <- paste('SELECT', fields, 'FROM', tbl)
    dbGetQuery(conn, sql)
}

setMethod("pmindex", "GenericPDInfo",
          function(object, subset=NULL, target=NULL) {
              ## target:
              ## - mps0: probes (without duplicates / pmfeature)
              ## - mpsK: K>= 1 (pmindex at the mpsK summary-level)
              if( substr(target, 1, 3) != 'mps')
                  stop("'target' for GenericArray-types must be in the 'mpsK'-format")
              if (target=='mps0'){
                  pmi <- getFeatureInfo(object, 'fid', type='pm')[[1]]
                  if (!is.null(subset))
                      warning("Subset not available for 'mps0'. Returning everything.")
              } else {
                  info <- getMPSInfo(object, substr(target, 4, 4), 'fid', type='pm')
                  ## info will have: fsetid, man_fsetid, fid
                  if (!is.null(subset))
                      info <- subset(info, man_fsetid %in% subset)
                  pmi <- info$fid
                  names(pmi) <- info$man_fsetid
                  rm(info)
              }
              sort(pmi)
          })


setMethod('pmindex', 'GenericFeatureSet',
          function(object, subset=NULL, target=NULL) {
              pmindex(get(annotation(object)), subset, target)
          })

setMethod("pm", "GenericFeatureSet",
          function(object, subset=NULL, target='mps1'){
            theClass <- class(exprs(object))
            pmi <- pmindex(object, subset=subset, target=target)
            if ("matrix" %in% theClass){
              out <- exprs(object)[pmi,, drop=FALSE]
            }else if ("ff_matrix" %in% theClass){
              out <- ffSubset(rows=pmi, object=exprs(object),
                              prefix="pm-")
            }
            return(out)
          })

setReplaceMethod("pm", signature(object="GenericFeatureSet", subset='ANY', target='ANY', value="matrix"),
                 function(object, subset=NULL, target='mps1', value){
                   tmp <- exprs(object)
                   tmp[pmindex(object, subset=subset, target=target),] <- value
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setReplaceMethod("pm", signature(object="GenericFeatureSet", subset='ANY', target='ANY', value="ff_matrix"),
                 function(object, subset=NULL, target='mps1', value){
                   tmp <- exprs(object)
                   open(tmp)
                   open(value)
                   finalizer(value) <- "delete"
                   nc <- ncol(tmp)
                   pmi <- pmindex(object, subset=subset, target=target)
                   for (i in 1:nc)
                       tmp[pmi, i] <- value[,i]
                   close(value)
                   close(tmp)
                   rm(value)
                   assayDataElementReplace(object, "exprs", tmp)
                 })

setMethod("rma", "GenericFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, target='mps1'){
              ## bgcorrect and normalization must be at the non-duplicated probe level
              pmi <- pmindex(object, subset=subset, target='mps0')
              pm0 <- exprs(object)[pmi,, drop=FALSE]
              if (background)
                  pm0 <- backgroundCorrect(pm0, method='rma')
              if (normalize)
                  pm0 <- normalize(pm0, method='quantile')
              rownames(pm0) <- as.character(pmi)
              targetInfo <- getMPSInfo(get(annotation(object)), substr(target, 4, 4), 'fid', type='pm')
              ##idx <- match(targetInfo$fid, pmi)
              ##pm0[idx,, drop=FALSE]
              exprs <- basicRMA(pm0[as.character(targetInfo$fid),,drop=FALSE],
                                pnVec=targetInfo$man_fsetid, normalize=FALSE,
                                background=FALSE, verbose=TRUE)
              colnames(exprs) <- sampleNames(object)
              out <- new("ExpressionSet")
              slot(out, "assayData") <- assayDataNew(exprs=exprs)
              slot(out, "phenoData") <- phenoData(object)
              slot(out, "featureData") <- basicAnnotatedDataFrame(exprs, byrow=TRUE)
              slot(out, "protocolData") <- protocolData(object)
              slot(out, "annotation") <- slot(object, "annotation")
              if (!validObject(out))
                  stop("Resulting object is invalid.")
              return(out)
          })

