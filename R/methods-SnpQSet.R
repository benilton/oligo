setMethod("initialize", "SnpQSet",
          function(.Object,
                   assayData = assayDataNew(senseThetaA=senseThetaA,
                     senseThetaB=senseThetaB,
                     antisenseThetaA=antisenseThetaA,
                     antisenseThetaB=antisenseThetaB),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
                   antisenseThetaA=new("matrix"),
                   antisenseThetaB=new("matrix"),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
                                  assayData = assayDataNew(
                                    senseThetaA=senseThetaA,
                                    senseThetaB=senseThetaB,
                                    antisenseThetaA=antisenseThetaA,
                                    antisenseThetaB=antisenseThetaB),
                                  phenoData=phenoData,
                                  experimentData=experimentData,
                                  annotation=annotation)
            .Object
          })

setValidity("SnpQSet",
            function(object)
            assayDataValidMembers(assayData(object),
                                  c("senseThetaA",
                                    "senseThetaB",
                                    "antisenseThetaA",
                                    "antisenseThetaB"))
            )


setMethod("senseThetaA", "SnpQSet", function(obj) assayData(obj)$senseThetaA)
setMethod("senseThetaB", "SnpQSet", function(obj) assayData(obj)$senseThetaB)
setMethod("antisenseThetaA", "SnpQSet", function(obj) assayData(obj)$antisenseThetaA)
setMethod("antisenseThetaB", "SnpQSet", function(obj) assayData(obj)$antisenseThetaB)
setMethod("getM", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)), ncol(antisenseThetaA(obj)), 2),
                         dimnames=list(rownames(antisenseThetaA(obj)), colnames(antisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- antisenseThetaA(obj)-antisenseThetaB(obj)
            tmp[,,2] <- senseThetaA(obj)-senseThetaB(obj)
            return(tmp)
          })
setMethod("getA", "SnpQSet",
          function(obj){
            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)), ncol(antisenseThetaA(obj)), 2),
                         dimnames=list(rownames(antisenseThetaA(obj)), colnames(antisenseThetaA(obj)),
                           c("antisense", "sense")))
            tmp[,,1] <- .5*(antisenseThetaA(obj)+antisenseThetaB(obj))
            tmp[,,2] <- .5*(senseThetaA(obj)+senseThetaB(obj))
            return(tmp)
          })
