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

getProbeIndex <- function(object, probeType='pm', target='core', subset, sortBy='fid'){
    probeType <- match.arg(probeType, getProbeTypes(object))
    target <- match.arg(target, getProbeTargets(object, probeType=probeType))
    sortBy <- match.arg(sortBy, c('fid', 'man_fsetid'))
    if (class(object) %in% c('ExonFeatureSet', 'GeneFeatureSet')){
        ## ST Arrays
        info <- switch(target,
               core=oligo:::getFidMetaProbesetCore(object, NULL),
               full=oligo:::getFidMetaProbesetFull(object, NULL),
               extended=oligo:::getFidMetaProbesetExtended(object, NULL),
               probeset=oligo:::getFidProbeset(object, NULL))
        names(info) <- c('fid', 'man_fsetid')
    }else{
        sql <- paste('SELECT fid, man_fsetid FROM',
                     paste(probeType, 'feature', sep=''))
        info <- dbGetQuery(db(object), sql)
    }
    sortCol <- which(names(info) == sortBy)
    if (!missing(subset))
        info <- subset(info, man_fsetid %in% subset)
    info[order(info[[sortCol]], info[[-sortCol]]), 'fid']
}

availProbeInfo <- function(object, probeType='pm', target='core'){
    ## FIXME: ST arrays have bg probes in pm tbl
    ## FIXME: ST arrays have targets in pm tbl (but same fields as pms?)
    dbListFields(db(object), paste(probeType, 'feature', sep=''))
}

getProbeInfo <- function(object, field, probeType='pm', target='core',
                         subset, sortBy='fid'){
    ## FIXME!!!
    probeType <- match.arg(probeType, getProbeTypes(object))
    target <- match.arg(target, getProbeTargets(object, probeType=probeType))
    avail2sort <- c('fid', 'man_fsetid')
    sortBy <- match.arg(sortBy, avail2sort)
    conn <- db(object)
    tbl <- match.arg(paste(probeType, 'feature', sep=''), dbListTables(conn))
    field <- match.arg(field, dbListFields(conn, tbl))
    if (class(object) %in% c('ExonFeatureSet', 'GeneFeatureSet')){
        ## FIXME: ST Arrays
        flds <- unique(paste(c('fid', 'fsetid as man_fsetid', field), collapse=', '))
        sql <- paste('SELECT', flds, 'FROM', tbl)
        info <- dbGetQuery(conn, sql)
    }else{
        flds <- unique(paste(c('fid', 'man_fsetid', field), collapse=', '))
        sql <- paste('SELECT', flds, 'FROM', tbl)
        info <- dbGetQuery(conn, sql)
    }
    if (!missing(subset))
        info <- subset(info, man_fsetid %in% subset)
    info[order(info[[sortBy]], info[[setdiff(avail2sort, sortBy)]]),]
}
