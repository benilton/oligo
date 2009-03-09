setMethod("rma", "ExonFeatureSet",
          function(object, background=TRUE, normalize=TRUE, subset=NULL){
            conn <- db(object)
            sql <- paste("SELECT fid, fsetid FROM pmfeature")
            featureInfo <- dbGetQuery(conn, sql)
            featureInfo <- featureInfo[order(featureInfo[["fsetid"]]),]
            rownames(featureInfo) <- NULL
            pms <- exprs(object)[featureInfo[["fid"]],, drop=FALSE]
            dimnames(pms) <- NULL
            theExprs <- basicRMA(pms,
                                 as.character(featureInfo[["fsetid"]]),
                                 length(unique(featureInfo[["fsetid"]])),
                                 normalize, background)
            rm(pms)

            ## Getting exon and gene info to add to featureData
            sql <- paste("SELECT fsetid, exon_id, transcript_cluster_id,",
                         "level, crosshyb_type, chrom",
                         "FROM featureSet")
            fromFSetTable <- dbGetQuery(conn, sql)
            sql <- paste("SELECT fsetid, group_concat(accession) as accessions",
                         "FROM fset2gene",
                         "INNER JOIN gene USING(gid)",
                         "GROUP BY fsetid")
            fromFSet2Gene <- dbGetQuery(conn, sql)
            fromFSet2Gene[["accessions"]] <- sapply(strsplit(fromFSet2Gene[["accessions"]], ","),
                                                    function(x) paste(sort(x), collapse=","))
            probesets <- merge(fromFSetTable, fromFSet2Gene, by="fsetid", all.x=TRUE)
            rm(sql, fromFSetTable, fromFSet2Gene)
            idx <- probesets[["fsetid"]] %in% unique(featureInfo[["fsetid"]])
            probesets <- probesets[idx,]
            rm(idx)
            
            descs <- c("probeset id", "exon id", "transcript cluster id",
                       "level", "cross hybridization type", "chromosome",
                       "accessions")
            vmd <- data.frame(labelDescription=descs)
            rm(descs)
            rownames(probesets) <- as.character(probesets[["fsetid"]])
            tmp <- new("AnnotatedDataFrame", data=probesets,
                       varMetadata=vmd)
            rm(vmd, probesets)

            out <- new("ExpressionSet",
                       featureData=tmp,
                       phenoData=phenoData(object),
                       annotation=annotation(object),
                       experimentData=experimentData(object),
                       exprs=theExprs)
            rm(tmp)
            return(out)
          })
