library(oligo)
xysPath <- system.file("extdata", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
summarized <- rma(ngsExpressionFeatureSet)
rm(ngsExpressionFeatureSet)
set.seed(1)
i0 <- sample(nrow(summarized), 1000)
save(i0, file='i0.rda')
rma_ref0 <- exprs(summarized)[i0,]
save(rma_ref0, file='rma_ref0.rda')

