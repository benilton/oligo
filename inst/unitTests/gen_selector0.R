pkgs <- c('pd.huex.1.0.st.v2', 'pd.hugene.1.0.st.v1')
getRandom <- function(dt){
    set.seed(1)
    i <- sort(sample(nrow(dt), 100))
    d0 <- dt[i,]
    rownames(d0) <- NULL
    list(idx=i, data=d0)
}

genFidCore <- function(pkg){
    library(pkg, character.only=TRUE)
    object <- get(pkg)
    conn <- db(object)
    sql <- "SELECT fid, meta_fsetid as man_fsetid FROM pmfeature INNER JOIN core_mps USING(fsetid) ORDER BY fid"
    featureInfo <- dbGetQuery(conn, sql)
    featureInfo$man_fsetid <- as.character(featureInfo$man_fsetid)
    featureInfo
}

genFidProbeset <- function(pkg){
    library(pkg, character.only=TRUE)
    object <- get(pkg)
    conn <- db(object)
    sql <- "SELECT fid, fsetid as man_fsetid FROM pmfeature ORDER BY fid"
    featureInfo <- dbGetQuery(conn, sql)
    featureInfo$man_fsetid <- as.character(featureInfo$man_fsetid)
    featureInfo
}

genAntigenomic <- function(pkg){
    library(pkg, character.only=TRUE)
    object <- get(pkg)
    conn <- db(object)
    sql <- 'SELECT type FROM type_dict WHERE type_id="control->bgp->antigenomic"'
    bgCode <- dbGetQuery(conn, sql)
    sql <- 'SELECT fid, fsetid as man_fsetid FROM pmfeature INNER JOIN featureSet USING(fsetid) WHERE type = :type ORDER BY fid'
    featureInfo <- dbGetPreparedQuery(conn, sql, bgCode)
    featureInfo$man_fsetid <- as.character(featureInfo$man_fsetid)
    featureInfo
}

core <- lapply(lapply(pkgs, genFidCore), getRandom)
pset <- lapply(lapply(pkgs, genFidProbeset), getRandom)
agen <- lapply(lapply(pkgs, genAntigenomic), getRandom)
core0 <- unlist(lapply(core, '[', 'data'), recursive=FALSE)
pset0 <- unlist(lapply(pset, '[', 'data'), recursive=FALSE)
agen0 <- unlist(lapply(agen, '[', 'data'), recursive=FALSE)
icore0 <- unlist(lapply(core, '[', 'idx'), recursive=FALSE)
ipset0 <- unlist(lapply(pset, '[', 'idx'), recursive=FALSE)
iagen0 <- unlist(lapply(agen, '[', 'idx'), recursive=FALSE)
names(core0) <- names(pset0) <- names(agen0) <- pkgs
names(icore0) <- names(ipset0) <- names(iagen0) <- pkgs
save(core0, pset0, agen0, icore0, ipset0, iagen0, file='fids_ref0.rda')
