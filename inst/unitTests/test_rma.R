rma_test <- function(){
  set.seed(1)
  xysPath <- system.file("extdata", package="maqcExpression4plex")
  xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
  ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
  summarized <- rma(ngsExpressionFeatureSet)
  rm(ngsExpressionFeatureSet)
  t0 <- exprs(summarized)[sample(nrow(summarized), 1000),]
  rm(summarized)
  load('rma_ref0.rda')
  all.equal(rma_ref0, t0)
}
