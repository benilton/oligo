## Boxplot
## begining of boxplot.R - by BC Th Jul 28, 2005

if( is.null(getGeneric("boxplot")))
  setGeneric("boxplot")

setMethod("boxplot",signature(x="oligoBatch"),
          function(x,which=c("both","pm","mm"),range=0,...){
            which <- match.arg(which,c("both","pm","mm"))
            tmp <- description(x)
            if (is(tmp, "MIAME")) main <- tmp@title

            tmp <- unlist(featureIndex(x,which))
            tmp <- tmp[seq(1,length(tmp),len=5000)]

            boxplot(data.frame(log2(exprs(x)[tmp,])),main=main,range=range, ...)
          })

## end of boxplot.R
