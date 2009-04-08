### For Tiling2

setMethod("pm", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- pmindex(object, subset=subset)
            pm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            pm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(pm1)
            out <- array(c(pm1, pm2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(pm1), colnames(pm1), c("channel1", "channel2"))
            return(out)
          })

setMethod("mm", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- mmindex(object, subset=subset)
            mm1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            mm2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(mm1)
            out <- array(c(mm1, mm2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(mm1), colnames(mm1), c("channel1", "channel2"))
            return(out)
          })

setMethod("bg", "TilingFeatureSet2",
          function(object, subset=NULL){
            idx <- bgindex(object, subset=subset)
            bg1 <- assayDataElement(object, "channel1")[idx,, drop=FALSE]
            bg2 <- assayDataElement(object, "channel2")[idx,, drop=FALSE]
            dims <- dim(bg1)
            out <- array(c(bg1, bg2), dim=c(dims, 2))
            dimnames(out) <- list(rownames(bg1), colnames(bg1), c("channel1", "channel2"))
            return(out)
          })


setReplaceMethod("pm", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- pmindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })

setReplaceMethod("mm", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- mmindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })

setReplaceMethod("bg", signature(object="TilingFeatureSet2", value="array"),
                 function(object, value){
                   idx <- bgindex(object)
                   tmp <- assayDataElement(object, "channel1")
                   tmp[idx,] <- value[,,1]
                   out <- assayDataElementReplace(object, "channel1", tmp)
                   tmp <- assayDataElement(out, "channel2")
                   tmp[idx,] <- value[,,2]
                   assayDataElementReplace(out, "channel2", tmp)
                 })
