setMethod("initialize", "DBPDInfo",
          function(.Object, ...) {
            .Object <- callNextMethod()
            tInfo <- dbGetQuery(db(.Object), "select * from table_info")
            .Object@tableInfo <- tInfo
            .Object
          })

setMethod("manufacturer", "PDInfo",
          function(object) object@manufacturer)

setMethod("genomeBuild", "PDInfo",
          function(object) object@genomebuild)

setMethod("db", signature(object="DBPDInfo"),
          function(object) object@getdb())

setMethod("geometry", "PDInfo",
          function(object) object@geometry)

setMethod("nrow", "platformDesign", function(x) x@nrow)
setMethod("ncol", "platformDesign", function(x) x@ncol)
## XXX: this should be type(object), but I ran into trouble
## because 'type' is somehow magic.  Didn't have time to
## track it down, so we're using kind for now.
setMethod("kind", "platformDesign", function(object) object@type)

setMethod("nProbes", "AffySNPPDInfo",
          function(object) {
            ## Note: does not include QC probes
            probeTables <- c("mmfeature", "pmfeature")
            sum(subset(object@tableInfo,
                       tbl %in% probeTables, row_count))
          })

setMethod("pmindex", "AffySNPPDInfo",
          function(object, subset=NULL) {
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

setMethod("kind", "AffyExpressionPDInfo",
          function(object) {
              "expression"
          })

setMethod("kind", "AffySNPPDInfo",
          function(object) {
              "SNP"
          })

setMethod("kind", "AffySNPCNVPDInfo",
          function(object) {
              "SNPCNV"
          })

setMethod("kind", "AffyGenePDInfo",
          function(object) {
              "Gene"
          })

setMethod("kind", "ExpressionPDInfo",
          function(object) {
            "expression"
          })

setMethod("kind", "TilingPDInfo",
          function(object) {
            "tiling"
          })

setMethod("probeNames", "AffySNPPDInfo",
          function(object, subset=NULL) {
            sql <- "select man_fsetid, fid from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            tmp[order(tmp$fid, tmp$man_fsetid), "man_fsetid"]
          })

## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("geneNames", "AffySNPPDInfo",
          function(object) {
              unique(probeNames(object))
          })

setMethod("pmSequence", "AffySNPPDInfo",
          function(object){
            sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
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

### For Expression

setMethod("nProbes", "ExpressionPDInfo",
          function(object) {
            ## Note: does not include QC probes
            probeTables <- c("mmfeature", "pmfeature")
            sum(subset(object@tableInfo,
                       tbl %in% probeTables, row_count))
          })

setMethod("pmindex", "ExpressionPDInfo",
          function(object, subset=NULL) {
            ## might improve by telling RSQLite how
            ## many rows we will fetch?  same for mmindex.
            dbGetQuery(db(object),
                       "select fid from pmfeature")[[1]]
          })

setMethod("mmindex", "ExpressionPDInfo",
          function(object) {
            dbGetQuery(db(object),
                       "select fid from mmfeature")[[1]]
          })

setMethod("probeNames", "ExpressionPDInfo",
          function(object, subset=NULL) {
            sql <- "select man_fsetid, fid from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            tmp[order(tmp$fid, tmp$man_fsetid), "man_fsetid"]
          })

## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("geneNames", "ExpressionPDInfo",
          function(object) {
              unique(probeNames(object))
          })

setMethod("pmSequence", "ExpressionPDInfo",
          function(object){
            sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmPosition", "ExpressionPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })


### For Tiling

setMethod("nProbes", "TilingPDInfo",
          function(object) {
            ## Note: does not include QC probes
            probeTables <- c("mmfeature", "pmfeature")
            sum(subset(object@tableInfo,
                       tbl %in% probeTables, row_count))
          })

setMethod("pmindex", "TilingPDInfo",
          function(object, subset=NULL) {
            ## might improve by telling RSQLite how
            ## many rows we will fetch?  same for mmindex.
            dbGetQuery(db(object),
                       "select fid from pmfeature")[[1]]
          })

setMethod("mmindex", "TilingPDInfo",
          function(object) {
            dbGetQuery(db(object),
                       "select fid from mmfeature")[[1]]
          })

setMethod("probeNames", "TilingPDInfo",
          function(object, subset=NULL) {
            sql <- "select man_fsetid, fid from featureSet, pmfeature where pmfeature.fsetid=featureSet.fsetid"
            tmp <- dbGetQuery(db(object), sql)
            tmp[order(tmp$fid, tmp$man_fsetid), "man_fsetid"]
          })

## FIXME: this method should be renamed!
## FIXME: this should query the featureSet table (much faster i think)
setMethod("geneNames", "TilingPDInfo",
          function(object) {
              unique(probeNames(object))
          })

setMethod("pmSequence", "TilingPDInfo",
          function(object){
            sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })

setMethod("pmPosition", "TilingPDInfo",
          function(object){
            sql <- "select position from pmfeature order by pmfeature.fid"
            dbGetQuery(db(object), sql)[[1]]
          })


### AffySNPCNV

setMethod("pmSequence", "AffySNPCNVPDInfo",
          function(object, probes.type="snp"){
            if (probes.type == "snp"){
              sql <- "select seq from sequence, pmfeature where pmfeature.fid=sequence.fid order by pmfeature.fid"
            }else if (probes.type == "cnv"){
              sql <- "SELECT seq FROM sequenceCNV, pmfeatureCNV WHERE pmfeatureCNV.fid=sequenceCNV.fid ORDER BY pmfeatureCNV.fid"
            }
            dbGetQuery(db(object), sql)[[1]]
          })
