### rearranged by function

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

setMethod("geometry", "DBPDInfo", ## changed signature from PDInfo to DBPDInfo
          function(object) object@geometry)
  
setMethod("db", signature(object="DBPDInfo"),
          function(object) object@getdb())

setMethod("nProbes", "DBPDInfo",
          function(object) {
            ## Note: does not include QC probes
            conn <- db(object)
            sql <- paste("SELECT row_count FROM table_info WHERE",
                         "tbl IN ('mmfeature', 'pmfeature', 'pmfeatureCNV')")
            sum(dbGetQuery(conn, sql)[[1]])
          })
  
## need to implement subset here, subset assumed to be Transcript Identifiers? MS
## PM
setMethod("pmindex", "AffySNPPDInfo",
          function(object, subset=NULL) { 
            ## might improve by telling RSQLite how
            ## many rows we will fetch?  same for mmindex.
            dbGetQuery(db(object),
                       "select fid from pmfeature")[[1]]
          })
## MM
setMethod("mmindex", "AffySNPPDInfo",
          function(object, subset=NULL) { ## has missing subset argument : MA
            dbGetQuery(db(object),
                       "select fid from mmfeature")[[1]]
          })

setMethod("kind", "AffySNPPDInfo",
		  function(object) {
			  "SNP"
		  })
  
setMethod("kind", "AffyExpressionPDInfo",
          function(object) {
              "expression"
          })


setMethod("kind", "AffySNPCNVPDInfo",
          function(object) {
              "SNPCNV"
          })
 
setMethod("kind", "AffyGenePDInfo",
          function(object) {
              "gene" ### changed to lower case to make consistant with others
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

## setMethod("nProbes", "ExpressionPDInfo",
##           function(object) {
##             ## Note: does not include QC probes
##             probeTables <- c("mmfeature", "pmfeature")
##             sum(subset(object@tableInfo,
##                        tbl %in% probeTables, row_count))
##           })

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

## setMethod("nProbes", "TilingPDInfo",
##           function(object) {
##             ## Note: does not include QC probes
##             probeTables <- c("mmfeature", "pmfeature")
##             sum(subset(object@tableInfo,
##                        tbl %in% probeTables, row_count))
##           })

setMethod("pmindex", "TilingPDInfo",
          function(object, subset=NULL) {
            ## might improve by telling RSQLite how
            ## many rows we will fetch?  same for mmindex.
            dbGetQuery(db(object),
                       "select fid from pmfeature")[[1]]
          })

setMethod("mmindex", "TilingPDInfo",
          function(object, subset=NULL) {
            if(!is.null(subset)) message("subset not yet implemented")
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
