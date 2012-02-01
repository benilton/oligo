test_rma <- function(){
  library(oligo)
  library(RUnit)
  xysPath <- system.file("extdata", package="maqcExpression4plex")
  xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
  ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
  summarized <- rma(ngsExpressionFeatureSet)
  rm(ngsExpressionFeatureSet)
##  set.seed(1)
##  idx <- sample(nrow(summarized), 1000)
load(system.file('unitTests', 'i0.rda', package='oligo'))
  t0 <- exprs(summarized)[i0,]
  rm(summarized)
  load(system.file('unitTests', 'rma_ref0.rda', package='oligo'))
  checkEquals(rma_ref0, t0)
}
