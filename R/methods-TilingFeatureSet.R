setMethod("pm", "TilingFeatureSet",
          function(object, subset=NULL){
            idx <- pmindex(object, subset=subset)
            if (channelNames(object) %in% c("channel1", "channel2")){
              pm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
              pm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
              dims <- dim(pm1)
              out <- array(c(pm1, pm2), dim=c(dims, 2))
              dimnames(out) <- list(rownames(pm1), colnames(pm1), c("channel1", "channel2"))
              return(out)
            }else{
              callNextMethod()
            }
          })

setMethod("mm", "TilingFeatureSet",
          function(object, subset=NULL){
            idx <- mmindex(object, subset=subset)
            if (channelNames(object) %in% c("channel1", "channel2")){
              mm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
              mm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
              dims <- dim(mm1)
              out <- array(c(mm1, mm2), dim=c(dims, 2))
              dimnames(out) <- list(rownames(mm1), colnames(mm1), c("channel1", "channel2"))
              return(out)
            }else{
              callNextMethod()
            }
          })

setMethod("bg", "TilingFeatureSet",
          function(object, subset=NULL){
            idx <- bgindex(object, subset=subset)
            if (channelNames(object) %in% c("channel1", "channel2")){
              bg1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
              bg2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
              dims <- dim(bg1)
              out <- array(c(bg1, bg2), dim=c(dims, 2))
              dimnames(out) <- list(rownames(bg1), colnames(bg1), c("channel1", "channel2"))
              return(out)
            }else{
              callNextMethod()
            }
          })


setReplaceMethod("pm", signature(object="TilingFeatureSet", value="array"),
                 function(object, value){
                   idx <- pmindex(object)
                   if (channelNames(object) %in% c("channel1", "channel2")){
                     tmp <- assayDataElement(object, "channel1")
                     tmp[idx,] <- value[,,1]
                     out <- assayDataElementReplace(object, "channel1", tmp)
                     tmp <- assayDataElement(out, "channel2")
                     tmp[idx,] <- value[,,2]
                     assayDataElementReplace(out, "channel2", tmp)
                   }else{
                     callNextMethod()
                   }
                 })

setReplaceMethod("mm", signature(object="TilingFeatureSet", value="array"),
                 function(object, value){
                   idx <- mmindex(object)
                   if (channelNames(object) %in% c("channel1", "channel2")){
                     tmp <- assayDataElement(object, "channel1")
                     tmp[idx,] <- value[,,1]
                     out <- assayDataElementReplace(object, "channel1", tmp)
                     tmp <- assayDataElement(out, "channel2")
                     tmp[idx,] <- value[,,2]
                     assayDataElementReplace(out, "channel2", tmp)
                   }else{
                     callNextMethod()
                   }
                 })

setReplaceMethod("bg", signature(object="TilingFeatureSet", value="array"),
                 function(object, value){
                   idx <- bgindex(object)
                   if (channelNames(object) %in% c("channel1", "channel2")){
                     tmp <- assayDataElement(object, "channel1")
                     tmp[idx,] <- value[,,1]
                     out <- assayDataElementReplace(object, "channel1", tmp)
                     tmp <- assayDataElement(out, "channel2")
                     tmp[idx,] <- value[,,2]
                     assayDataElementReplace(out, "channel2", tmp)
                   }else{
                     callNextMethod()
                   }
                 })

setMethod("getContainer", "TilingFeatureSet",
          function(object, probeType=c("pm", "bg")){
            probeType <- match.arg(probeType)
            tbl <- ifelse(probeType == "pm", "pmfeature", "bgfeature")
            conn <- db(object)
            sql <- paste("SELECT type FROM", tbl,
                         "INNER JOIN featureSet",
                         "USING (fsetid)",
                         "ORDER BY fid")
            dbGetQuery(conn, sql)[[1]]
          })

setMethod("getM", "TilingFeatureSet",
          function(object){
            chn <- channelNames(object)
            if (!(chn %in% c("channel1", "channel2")))
              stop("Object does not have 'channel1' and 'channel2'.")
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))
            lc1-lc2
          })

setMethod("getA", "TilingFeatureSet",
          function(object){
            chn <- channelNames(object)
            if (!(chn %in% c("channel1", "channel2")))
              stop("Object does not have 'channel1' and 'channel2'.")
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))
            (lc1+lc2)/2
          })

setMethod("pmPosition", "TilingFeatureSet",
          function(object){
            conn <- db(object)
            tmp <- dbGetQuery(conn, "SELECT fid, position FROM pmfeature")
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })
