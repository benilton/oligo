## Boxplot
## begining of boxplot.R - by BC Th Jul 28, 2005

## FIX ME!!!
## BC: it seems that the final code (after building)
##     is created using the .R files in ascending order.
##     boxplot() uses an oligoBatch argument and
##     to get rid of some warning messages I defined the
##     class in the first file that appears in the directory
setClass("oligoBatch",
         representation(manufacturer="character",
                        platform="character"),
         contains="eSet")

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
