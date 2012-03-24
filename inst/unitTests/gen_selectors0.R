pkgs <- c('pd.huex.1.0.st.v2', 'pd.hugene.1.0.st.v1')
getRandom <- function(dt){
    set.seed(1)
    i <- sort(sample(nrow(dt), 100))
    d0 <- dt[i,]
    rownames(d0) <- NULL
    d0
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
core0 <- lapply(lapply(pkgs, genFidCore), getRandom)
pset0 <- lapply(lapply(pkgs, genFidProbeset), getRandom)
agen0 <- lapply(lapply(pkgs, genAntigenomic), getRandom)
names(core0) <- names(pset0) <- names(agen0) <- pkgs
save(core0, pset0, agen0, file='fids_ref0.rda')
