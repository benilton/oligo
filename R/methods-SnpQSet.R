setMethod("senseThetaA", "SnpQSet", function(object) assayData(object)$senseThetaA)
setMethod("senseThetaB", "SnpQSet", function(object) assayData(object)$senseThetaB)
setMethod("antisenseThetaA", "SnpQSet", function(object) assayData(object)$antisenseThetaA)
setMethod("antisenseThetaB", "SnpQSet", function(object) assayData(object)$antisenseThetaB)

setMethod("getM", "SnpQSet",
          function(object){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(object)),
                               ncol(antisenseThetaA(object)), 2),
                         dimnames=list(rownames(antisenseThetaA(object)),
                           colnames(antisenseThetaA(object)),
                           c("antisense", "sense")))
            tmp[,,1] <- antisenseThetaA(object)-antisenseThetaB(object)
            tmp[,,2] <- senseThetaA(object)-senseThetaB(object)
            return(tmp)
          })

setMethod("getA", "SnpQSet",
          function(object){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(object)),
                               ncol(antisenseThetaA(object)), 2),
                         dimnames=list(rownames(antisenseThetaA(object)),
                           colnames(antisenseThetaA(object)),
                           c("antisense", "sense")))
            tmp[,,1] <- .5*(antisenseThetaA(object)+antisenseThetaB(object))
            tmp[,,2] <- .5*(senseThetaA(object)+senseThetaB(object))
            return(tmp)
          })

setMethod("db", "SnpQSet", function(object) db(get(annotation(object))))

setMethod("plotM", c("SnpQSet", "integer"),
          function(object, i, ...){
            mm <- featureNames(object[i,])
            plot(getM(object[i,])[,,], main=mm, ...)})

setMethod("plotM", c("SnpQSet", "numeric"),
          function(object, i, ...)
          plotM(object, as.integer(i), ...))

setMethod("plotM", c("SnpQSet", "character"),
          function(object, i, ...){
            ii <- as.integer(which(featureNames(object) == i))
            plotM(object, ii, ...)
          })

##

setMethod("thetaA", "SnpCnvQSet", function(object) assayDataElement(object, "thetaA"))
setMethod("thetaB", "SnpCnvQSet", function(object) assayDataElement(object, "thetaB"))
setMethod("getM", "SnpCnvQSet", function(object) thetaA(object)-thetaB(object))
setMethod("getA", "SnpCnvQSet", function(object) (.5*thetaA(object)+.5*thetaB(object)))
setMethod("db", "SnpCnvQSet", function(object) db(get(annotation(object))))
