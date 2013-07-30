setClass('oligoPLM',
         representation(chip.coefs='matrix',
                        probe.coefs='numeric',
                        weights='matrix',
                        residuals='matrix',
                        se.chip.coefs='matrix',
                        se.probe.coefs='numeric',
                        residualSE='numeric',
                        geometry='numeric',
                        method='character',
                        manufacturer='character',
                        annotation='character',
                        narrays='integer',
                        nprobes='integer',
                        nprobesets='integer')
         )

setMethod('ncol', 'oligoPLM',
          function(object){
              object@narrays
          })

setGeneric('nprobes', function(object) standardGeneric('nprobes'))
setMethod('nprobes', 'oligoPLM',
          function(object){
              object@nprobes
          })

setGeneric('nprobesets', function(object) standardGeneric('nprobesets'))
setMethod('nprobesets', 'oligoPLM',
          function(object){
              object@nprobesets
          })

setGeneric('coefs', function(object) standardGeneric('coefs'))
setMethod('coefs', 'oligoPLM',
          function(object){
              object@chip.coefs
          })

setGeneric('coefs<-', function(object, value) standardGeneric('coefs<-'))
setReplaceMethod('coefs', 'oligoPLM',
                 function(object, value){
                     object@chip.coefs <- value
                     object
                 })

setGeneric('coefs.probe', function(object) standardGeneric('coefs.probe'))
setMethod('coefs.probe', 'oligoPLM',
          function(object){
              object@probe.coefs
          })

setGeneric('coefs.probe<-', function(object, value) standardGeneric('coefs.probe<-'))
setReplaceMethod('coefs.probe', 'oligoPLM',
                 function(object, value){
                     object@probe.coefs <- value
                     object
                 })

##setGeneric('weights', function(object) standardGeneric('weights'))
setMethod('weights', 'oligoPLM',
          function(object, ...){
              object@weights
          })

setGeneric('weights<-', function(object, value) standardGeneric('weights<-'))
setReplaceMethod('weights', 'oligoPLM',
                 function(object, value){
                     object@weights <- value
                     object
                 })

setGeneric('resids', function(object) standardGeneric('resids'))
setMethod('resids', 'oligoPLM',
          function(object){
              object@residuals
          })

setGeneric('resids<-', function(object, value) standardGeneric('resids<-'))
setReplaceMethod('resids', 'oligoPLM',
                 function(object, value){
                     object@residuals <- value
                     object
                 })

setGeneric('se', function(object) standardGeneric('se'))
setMethod('se', 'oligoPLM',
          function(object){
              object@se.chip.coefs
          })

setGeneric('se<-', function(object, value) standardGeneric('se<-'))
setReplaceMethod('se', 'oligoPLM',
                 function(object, value){
                     object@se.chip.coefs <- value
                     object
                 })

setGeneric('se.probe', function(object) standardGeneric('se.probe'))
setMethod('se.probe', 'oligoPLM',
          function(object){
              object@se.probe.coefs
          })

setGeneric('se.probe<-', function(object, value) standardGeneric('se.probe<-'))
setReplaceMethod('se.probe', 'oligoPLM',
                 function(object, value){
                     object@se.probe.coefs <- value
                     object
                 })

setGeneric('residualSE', function(object) standardGeneric('residualSE'))
setMethod('residualSE', 'oligoPLM',
          function(object){
              object@residualSE
          })

setGeneric('residualSE<-', function(object, value) standardGeneric('residualSE<-'))
setReplaceMethod('residualSE', 'oligoPLM',
                 function(object, value){
                     object@residualSE <- value
                     object
                 })

setGeneric('geometry', function(object) standardGeneric('geometry'))
setMethod('geometry', 'oligoPLM',
          function(object){
              object@geometry
          })

setGeneric('geometry<-', function(object, value) standardGeneric('geometry<-'))
setReplaceMethod('geometry', 'oligoPLM',
                 function(object, value){
                     object@geometry <- value
                     object
                 })

setGeneric('method', function(object) standardGeneric('method'))
setMethod('method', 'oligoPLM',
          function(object){
              object@method
          })

setGeneric('method<-', function(object, value) standardGeneric('method<-'))
setReplaceMethod('method', 'oligoPLM',
                 function(object, value){
                     object@method <- value
                     object
                 })

setGeneric('manufacturer', function(object) standardGeneric('manufacturer'))
setMethod('manufacturer', 'oligoPLM',
          function(object){
              object@manufacturer
          })

setGeneric('manufacturer<-', function(object, value) standardGeneric('manufacturer<-'))
setReplaceMethod('manufacturer', 'oligoPLM',
                 function(object, value){
                     object@manufacturer <- value
                     object
                 })

setGeneric('annotation', function(object) standardGeneric('annotation'))
setMethod('annotation', 'oligoPLM',
          function(object){
              object@annotation
          })

setGeneric('annotation<-', function(object, value) standardGeneric('annotation<-'))
setReplaceMethod('annotation', 'oligoPLM',
                 function(object, value){
                     object@annotation <- value
                     object
                 })

setMethod("show", "oligoPLM",
          function(object) {

              ## Add preprocessing info
              message("Probe Level Model")
              message("Method..............: ", object@method)
              message("# Samples...........: ", ncol(object))
              message("# Probes (processed): ", nprobes(object))
              message("# Probesets.........: ", nprobesets(object))
              message("Annotation..........: ", annotation(object))
          })

## fix names(theMat)
setMethod("boxplot",signature(x="oligoPLM"),
          function(x, type=c("NUSE", "RLE", "weights","resids"), col=darkColors(ncol(x)), range=0, ylim, ...){
            type <- match.arg(type)
            if (type == 'NUSE'){
                theMat <- NUSE(x, type='values')
                candYL <- c(.95, 1.10)
            }else if (type == 'RLE'){
                theMat <- RLE(x, type='values')
                candYL <- c(-.75, .75)
            }else if (type == 'weights'){
                theMat <- weights(x)
                theMat <- theMat[!is.na(theMat[,1]),,drop=FALSE]
                candYL <- c(0, 1)
            }else{
                theMat <- resids(x)
                theMat <- theMat[!is.na(theMat[,1]),,drop=FALSE]
                candYL <- c(-1, 1)
            }
            if (missing(ylim)) ylim <- candYL
            theMat <- as.data.frame(theMat)
            boxplot(theMat, col=col, range=range, ylim=ylim, ...)
          })

