test_rma <- function(){
  library(oligo)
  library(RUnit)
  message('Getting sample dataset')
  xysPath <- system.file("extdata", package="maqcExpression4plex")
  xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
  ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
  message('Running RMA')
  summarized <- rma(ngsExpressionFeatureSet)
  rm(ngsExpressionFeatureSet)
  load(system.file('unitTests', 'i0.rda', package='oligo'))
  t0 <- exprs(summarized)[i0,]
  rm(summarized)
  message('Getting reference results')
  load(system.file('unitTests', 'rma_ref0.rda', package='oligo'))
  checkEquals(rma_ref0, t0)
}
