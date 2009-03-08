setMethod("rma", "ExonFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL, level="exon", verbose=TRUE){
            level <- match.arg(level, c("exon", "gene"))
            conn <- db(object)

            if (verbose) message("RMA on Exon arrays at the ", level, " level.")

            ## The grouping variable (exon/gene) below will be called 'grp'
            ## Select featuresets with the following conditions:
            ##  level = core (1)  (from featureSet)
            ##  xhyb = unique (1) (from featureSet)
            ##  probeset maps to only one gene accession
            if (level == "exon"){
              sql <- paste("SELECT fid, exon_id as grp",
                           "FROM pmfeature",
                           "INNER JOIN featureSet",
                           "USING(fsetid)",
                           "WHERE level=1 AND crosshyb_type=1")
              featureInfo <- dbGetQuery(conn, sql)
              rm(sql)
              featureInfo <- featureInfo[order(featureInfo[["grp"]]),]
              rownames(featureInfo) <- NULL
              featureInfo[["grp"]] <- as.character(featureInfo[["grp"]])
            } else if (level == "gene"){
              sql <- paste("SELECT DISTINCT featureSet.fsetid, accession as grp",
                           "FROM featureSet, fset2gene, gene",
                           "WHERE level=1 AND crosshyb_type=1 AND",
                           "featureSet.fsetid=fset2gene.fsetid AND",
                           "fset2gene.gid=gene.gid",
                           "GROUP BY featureSet.fsetid",
                           "HAVING COUNT(DISTINCT accession)=1")
              coreGenesUnique <- dbGetQuery(conn, sql)
              coreGenesUnique <- coreGenesUnique[order(coreGenesUnique[["fsetid"]]),]
              rm(sql)

              pmProbesSql <- paste("SELECT fid, pmfeature.fsetid",
                                   "FROM pmfeature",
                                   "INNER JOIN featureSet",
                                   "USING(fsetid)",
                                   "WHERE level=1 AND crosshyb_type=1")
              pmProbes <- dbGetQuery(conn, pmProbesSql)
              rm(pmProbesSql)
              idx <- which(pmProbes[["fsetid"]] %in% coreGenesUnique[["fsetid"]])
              pmProbes <- pmProbes[idx,]
              rm(idx)
              featureInfo <- merge(pmProbes, coreGenesUnique,
                                   by.x="fsetid", by.y="fsetid")[, c("fid", "grp")]
              rm(coreGenesUnique, pmProbes)
            } else {
              stop("Unknown level selected: ", level)
            }

            ## Remove groups with less than 4 probes
            toRemove <- table(featureInfo[["grp"]])
            toRemove <- names(toRemove[toRemove < 4])
            if (verbose) message("Removing ", length(toRemove), " ", level, "s for not having at least 4 probes.")
            toRemove <- which(featureInfo[["grp"]] %in% toRemove)
            featureInfo <- featureInfo[-toRemove,]
            rm(toRemove)
            
            pnVec <- featureInfo[["grp"]]
            ngenes <- length(unique(pnVec))
            pms <- exprs(object)[featureInfo[["fid"]],, drop=FALSE]

            bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}

            exprs <- basicRMA(pms, pnVec, ngenes, normalize, background)
            
            out <- new("ExpressionSet",
                       phenoData=phenoData(object),
                       annotation=annotation(object),
                       experimentData=experimentData(object),
                       exprs=exprs)

            tmp <- preproc(object)
            tmp[["analysis"]] <- paste("RMA:", level, "level")
            preproc(out) <- tmp
            rm(tmp)
            return(out)
          })
