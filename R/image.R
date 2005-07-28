## image.R
## Thu, Jul 28, 2005 by BC

if( is.null(getGeneric("image")))
  setGeneric("image")

setMethod("image",signature(x="oligoBatch"),
          function(x, transfo=log, col=gray(c(0:64)/64),xlab="",ylab="", ...){
            scn <- prod(par("mfrow"))
            ask <- dev.interactive()
            which.plot <- 0

            ## x.pos <- (1:nrow(x)) - (1 + getOption("BioC")$affy$xy.offset)
            ## y.pos <- (1:ncol(x)) - (1 + getOption("BioC")$affy$xy.offset)

            x.pos <- (1:nrow(x)) - 1
            y.pos <- (1:ncol(x)) - 1

            for(i in 1:length(sampleNames(x))){
              which.plot <- which.plot+1;
              if(trunc((which.plot-1)/scn)==(which.plot-1)/scn && which.plot>1 && ask)  par(ask=TRUE)
              m <- exprs(x)[,i]
              if (is.function(transfo)) {
                m <- transfo(m)
              }
              m <- as.matrix(rev(as.data.frame(matrix(m, nrow=length(x.pos), ncol=length(y.pos)))))
              image(x.pos, y.pos, m,
                    col=col, main=sampleNames(x)[i],
                    xlab=xlab, ylab=ylab,,xaxt='n',
                      yaxt='n', ...)
              par(ask=FALSE)
            }
          })

