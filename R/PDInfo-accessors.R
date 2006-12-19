setMethod("initialize", "DBPDInfo",
          function(.Object) {
              .Object <- callNextMethod()
              tInfo <- dbGetQuery(db(.Object), "select * from table_info")
              .Object@tableInfo <- tInfo
              .Object
          })



setMethod("manufacturer", "PDInfo",
          function(object) object@manufacturer)

setMethod("genomeBuild", "PDInfo",
                    function(object) object@genomebuild)

setMethod("db", "DBPDInfo",
          function(object) object@getdb())


setMethod("nrow", "platformDesign", function(x) x@nrow)
setMethod("ncol", "platformDesign", function(x) x@ncol)
## XXX: this should be type(object), but I ran into trouble
## because 'type' is somehow magic.  Didn't have time to
## track it down, so we're using kind for now.
setMethod("kind", "platformDesign", function(object) object@type)


setMethod("listFeatureFields", "AffySNPPDInfo",
           function(object) {
               PM_TABLE <- "pmfeature"
               dbListFields(db(object), PM_TABLE)
           })

setMethod("listFeatureSetFields", "AffySNPPDInfo",
           function(object) {
               FSET_TABLE <- "featureSet"
               dbListFields(db(object), PM_TABLE)
           })

setMethod("nProbes", "AffySNPPDInfo",
          function(object) {
              ## Note: does not include QC probes
              probeTables <- c("mmfeature", "pmfeature")
              sum(subset(object@tableInfo,
                         tbl %in% probeTables, row_count))
          })

setMethod("pmindex", "AffySNPPDInfo",
          function(object) {
              ## might improve by telling RSQLite how
              ## many rows we will fetch?  same for mmindex.
              dbGetQuery(db(object),
                         "select fid from pmfeature")[[1]]
          })

setMethod("mmindex", "AffySNPPDInfo",
          function(object) {
              dbGetQuery(db(object),
                         "select fid from mmfeature")[[1]]
          })

allPMAllele <- function(map) {
    sqliteQuickColumn(db(map), "pmfeature", "allele")
}


allPMStrand <- function(map) {
    sqliteQuickColumn(db(map), "pmfeature", "strand")
}

pmIdsByAllele1 <- function(map, allele=c("A", "B")) {
    if (allele == "A")
      allPMIds(map)[!allPMAllele(map)]
    else
      allPMIds(map)[allPMAllele(map)]
}

pmIdsByAllele2 <- function(map, allele=c("A", "B")) {
    sql <- "select fid from pmfeature where allele = "
    if (allele == "A")
      sql <- paste(sql, "'1'")
    else
      sql <- paste(sql, "'0'")
    ans <- dbGetQuery(db(map), sql)[[1]]
    ans
}

pmFeatures <- function(map, featureSetIds) {
    fsetIds <- paste("'", featureSetIds, "'", sep="", collapse=",")
    sql <- paste("select man_fsetid, pmfeature.fid, pmfeature.strand,
pmfeature.allele, pmfeature.x, pmfeature.y from
pmfeature, featureSet where man_fsetid in (",
                 fsetIds, ") and featureSet.fsetid = pmfeature.fsetid")
    dbGetQuery(db(map), sql)
}

allFeatureSetIds <- function(map) {
    sqliteQuickColumn(db(map), "featureSet", "man_fsetid")
}


setMethod("kind", "AffySNPPDInfo",
          function(object) {
              "SNP"
          })

setMethod("featureSetNames", "AffySNPPDInfo",
          function(object) {
              ## FIXME, we may need to remove QC featureSets?
              sql <- "select man_fsetidid from featureSet"
              dbGetQuery(db(object), sql)[[1]]
          })
