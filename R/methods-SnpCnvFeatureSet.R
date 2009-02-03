## setMethod("pmindex", "AffySNPCNVPDInfo",
##           function(object, probes.type="snp") {
##             ## might improve by telling RSQLite how
##             ## many rows we will fetch?  same for mmindex.
##             if (probes.type == "snp"){
##               dbGetQuery(db(object),
##                          "select fid from pmfeature")[[1]]
##             }else if (probes.type == "cnv"){
##               dbGetQuery(db(object),
##                          "SELECT fid FROM pmfeatureCNV ORDER BY fid")[[1]]
##             }
##           })
## 
## setMethod("pmindex", "SnpCnvFeatureSet",
##           function(object, probes.type="snp"){
##             pmindex(getPlatformDesign(object), probes.type=probes.type)
##           })
## 
## setMethod("pm", "SnpCnvFeatureSet",
##           function(object, subset=NULL, probes.type="snp"){
##             if (!is.null(subset)) message("subset ignored (not implemented yet)")
##             tmp <- subBufferedMatrix(exprs(object), pmindex(object, probes.type=probes.type))
##             RowMode(tmp)
##             set.buffer.dim(tmp, 50000, 1)
##             return(tmp)
##           })
