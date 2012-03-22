######################################################
######################################################
### Tools for selecting probes and simplying tasks
### within oligo
### By: Benilton Carvalho - Jan/12
######################################################
######################################################

## FIXME: add a getProbeLevels??? (ST: bg/antigenomic/genomic/pm/etc)

getProbeTypes <- function(object){
    ## for exon/gene ST arrays, the BG probes are stored in the
    ## pmfeature as well... how to fix?
    substr(grep("^[a-z]{2}feature$", dbListTables(db(object)),
                value=TRUE), 1, 2)
}

getProbeTargets <- function(object, probeType='pm'){
    ## TODO: fix me!
    switch(probeType,
           pm=c('core', 'full', 'extended', 'probeset'),
           bg=c('genomic', 'antigenomic'))
}

## getProbeIndex <- function(object, probeType='pm', target='core', subset, sortBy='fid'){
##     probeType <- match.arg(probeType, getProbeTypes(object))
##     target <- match.arg(target, getProbeTargets(object, probeType=probeType))
##     sortBy <- match.arg(sortBy, c('fid', 'man_fsetid'))
##     if (class(object) %in% c('ExonFeatureSet', 'GeneFeatureSet')){
##         ## ST Arrays
##         info <- switch(target,
##                core=oligo:::getFidMetaProbesetCore(object, NULL),
##                full=oligo:::getFidMetaProbesetFull(object, NULL),
##                extended=oligo:::getFidMetaProbesetExtended(object, NULL),
##                probeset=oligo:::getFidProbeset(object, NULL))
##         names(info) <- c('fid', 'man_fsetid')
##     }else{
##         sql <- paste('SELECT fid, man_fsetid FROM',
##                      paste(probeType, 'feature', sep=''))
##         info <- dbGetQuery(db(object), sql)
##     }
##     sortCol <- which(names(info) == sortBy)
##     if (!missing(subset))
##         info <- subset(info, man_fsetid %in% subset)
##     info[order(info[[sortCol]], info[[-sortCol]]), 'fid']
## }

availProbeInfo <- function(object, probeType='pm', target='core'){
    if (!require(annotation(object), character.only=TRUE))
        stop("The annotation package '", annotation(object), "' is not available.")
    ## FIXME: ST arrays have bg probes in pm tbl
    ## FIXME: ST arrays have targets in pm tbl (but same fields as pms?)
    isST <- class(object) %in% c('ExonFeatureSet', 'GeneFeatureSet')
    conn <- db(object)
    probeTable <- paste(probeType, 'feature', sep='')
    ptFields <- dbListFields(conn, probeTable)
    probesetTable <- 'featureSet'
    psFields <- dbListFields(conn, probesetTable)
    if (isST & target!='probeset') psFields <- setdiff(psFields, 'transcript_cluster_id')
    out <- list(ptFields, psFields)
    names(out) <- c(probeTable, probesetTable)
    if (isST & target!='probeset'){
        mpsTable <- paste(target, 'mps', sep='_')
        mpFields <- dbListFields(conn, mpsTable)
        out[[3]] <- mpFields
        names(out) <- c(probeTable, probesetTable, mpsTable)
    }
    out
}

getProbeInfo <- function(object, field, probeType='pm', target='core',
                         subset, sortBy=c('fid', 'man_fsetid', 'none')){
    if (!require(annotation(object), character.only=TRUE))
        stop("The annotation package '", annotation(object), "' is not available.")
    sortBy <- match.arg(sortBy)
    conn <- db(object)

    ## With ST arrays:
    ## 1) fsetid is both fsetid and man_fsetid
    ## 2) transcript_cluster_id is in both *mps and featureSet tables
    ## 3) chrom/level/type_dict tables exist
    isST <- class(object) %in% c('ExonFeatureSet', 'GeneFeatureSet')

    if (missing(field)) field <- 'fid'

    probeTable <- paste(probeType, 'feature', sep='')
    if (isST & target!='probeset'){
        fields <- unique(c('fid', 'meta_fsetid as man_fsetid', field))
        fields[fields == 'fsetid'] <- 'pmfeature.fsetid'
        fields <- paste(fields, collapse=', ')
        mpsTable <- paste(target, 'mps', sep='_')
        fields <- gsub('transcript_cluster_id',
                       paste(mpsTable, '.transcript_cluster_id as transcript_cluster_id', sep=''),
                       fields)
        tables <- paste('pmfeature, featureSet,', mpsTable)
        sql <- paste('SELECT', fields, 'FROM', tables,
                     'WHERE pmfeature.fsetid=featureSet.fsetid AND',
                     paste('featureSet.fsetid=', mpsTable, '.fsetid', sep=''))
        rm(fields, mpsTable, tables)
    }else{
        fields <- unique(c('fid', 'man_fsetid', field))
        fields[fields == 'fsetid'] <- paste(probeTable, 'fsetid', sep='.')
        if (isST | class(object) == 'TilingFeatureSet')
            fields[fields == 'man_fsetid'] <- paste(probeTable, 'fsetid as man_fsetid', sep='.')
        fields <- paste(fields, collapse=', ')
        sql <- paste('SELECT', fields, 'FROM',
                     probeTable, 'INNER JOIN featureSet',
                     'USING(fsetid)')
        rm(fields)
    }

    info <- dbGetQuery(conn, sql)
    ## TODO: Add subset here!

    ## Getting data from dictionaries
    ## .. chrom: chromosome: mapping to proper chr ids
    ## .. level: ST arrays: bg/genomic/antigenomic
    ## .. type: ST arrays: core, extended, full
    field2dict <- c('chrom', 'level', 'type')
    addMerge <- field2dict[field2dict %in% field]
    if (length(addMerge) > 0)
        for (item in addMerge){
            dict <- paste(item, 'dict', sep='_')
            sql <- paste("SELECT * FROM",dict)
            info <- merge(info, dbGetQuery(conn, sql), all.x=TRUE, sort=FALSE)
            info[[item]] <- NULL
            rm(dict, sql)
        }

    ## This is to rename the dictionary columns that are
    ##   like chrom_id, type_id, level_id
    nmsOrig <- names(info)
    iTID <- grep('transcript_cluster_id', nmsOrig)
    origValue <- NULL
    if (length(iTID) > 0) origValue <- nmsOrig[iTID]
    nmsOrig <- gsub('\\_id$', '', nmsOrig)
    if (length(iTID) > 0) nmsOrig[iTID] <- origValue
    names(info) <- nmsOrig
    rm(nmsOrig, iTID, origValue)

    ## Sorting
    if (sortBy!='none'){
        i2 <- setdiff(c('fid', 'man_fsetid'), sortBy)
        info <- info[order(info[[sortBy]], info[[i2]]),]
        rownames(info) <- NULL
        rm(i2)
    }

    ## Make sure man_fsetid is character
    if ('man_fsetid' %in% names(info))
        if (!is.character(info$man_fsetid))
            info$man_fsetid <- as.character(info$man_fsetid)

    ## Return final object
    info
}
