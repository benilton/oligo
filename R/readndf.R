# Reading NimbleGen Data
# Author: Benilton Carvalho
# Date: Mar / 2005

readndf <- function(ndffile){
  header <- scan(ndffile,
                 what=character(0),
                 sep="\t",
                 nlines=1,
                 quiet=T)
  
  types <- list(PROBE_DESIGN_ID="character",
                CONTAINER="character",
                DESIGN_NOTE="character",
                SELECTION_CRITERIA="character",
                SEQ_ID="character",
                PROBE_SEQUENCE="character",
                MISMATCH="numeric",
                MATCH_INDEX="numeric",
                FEATURE_ID="numeric",
                ROW_NUM="numeric",
                COL_NUM="numeric",
                PROBE_CLASS="character",
                PROBE_ID="character",
                POSITION="numeric",
                DESIGN_ID="numeric",
                X="numeric",
                Y="numeric")
  whatToRead <- types[match(header,names(types))]
  
  # Reading the Nimblegen Description File
  # Using scan to make it faster.
  ndfdata <- scan(ndffile,what=whatToRead,skip=1,sep="\t",quiet=T)
  
  ndfdata$MISMATCH <- as.integer(ndfdata$MISMATCH)
  ndfdata$MATCH_INDEX <- as.integer(ndfdata$MATCH_INDEX)
  ndfdata$FEATURE_ID <- as.integer(ndfdata$FEATURE_ID)
  #ndfdata$ROW_NUM <- as.integer(ndfdata$ROW_NUM)
  #ndfdata$COL_NUM <- as.integer(ndfdata$COL_NUM)
  #ndfdata$POSITION <- as.integer(ndfdata$POSITION)
  ndfdata$DESIGN_ID <- as.integer(ndfdata$DESIGN_ID)
  
  ndfdata$X <- as.integer(ndfdata$X)
  ndfdata$Y <- as.integer(ndfdata$Y)
#  ndfdata$INDEX <- max(ndfdata$Y)*ndfdata$X + ndfdata$Y
  ndfdata$PROBESET <- as.factor(ndfdata$PROBE_ID)
#  ndfdata$PROBETYPE[ndfdata$MISMATCH == 0 & ndfdata$MATCH_INDEX!=0]<- "PM"
#  ndfdata$PROBETYPE[ndfdata$MISMATCH == 1 & ndfdata$MATCH_INDEX!=0] <- "MM"
#  ndfdata$PROBETYPE[ndfdata$PROBE_CLASS == "control"] <- "control"
  ndfdata$PROBETYPE[ndfdata$MISMATCH == 0 & ndfdata$PROBE_CLASS == "experimental"]<- "PM"
  ndfdata$PROBETYPE[ndfdata$MISMATCH >= 1 & ndfdata$PROBE_CLASS == "experimental"] <- "MM"
  ndfdata$PROBETYPE[ndfdata$PROBE_CLASS != "experimental"] <- "control"
  ndfdata$PROBETYPE <- as.factor(ndfdata$PROBETYPE)
  ndfdata$SEQUENCE <- ndfdata$PROBE_SEQUENCE
  ndfdata$FEATURE <- ndfdata$FEATURE_ID
  ndfdata$DESIGN <- factor(ndfdata$DESIGN_ID)
  ndfdata$PROBE_ID <- factor(ndfdata$PROBE_ID)
#  ndfdata$DESIGN_ID <- NULL
#  ndfdata$SEQ_ID <- NULL
#  ndfdata$FEATURE_ID <- NULL
#  ndfdata$PROBE_SEQUENCE <- NULL
#  ndfdata$MISMATCH <- NULL
#  ndfdata$PROBE_CLASS <- NULL
#  ndfdata$MATCH_INDEX <- NULL
#  ndfdata$COL_NUM <- NULL
#  ndfdata$ROW_NUM <- NULL
#  ndfdata$POSITION <- NULL
#  ndfdata$PROBE_ID <- NULL
#  ndfdata$DESIGN_NOT <- NULL
#  ndfdata$SELECTION_CRITERIA <- NULL
#  ndfdata$PROBE_DESIGN_ID <- NULL
#  ndfdata$DESIGN_NOTE <- NULL

  ndfdata2 <- data.frame(feature_id = ndfdata$FEATURE_ID,
                         X = ndfdata$X,
                         Y = ndfdata$Y,
                         feature_names = ndfdata$PROBE_ID,
                         sequence = ndfdata$PROBE_SEQUENCE,
                         feature_type = ndfdata$PROBETYPE,
                         mismatch = ndfdata$MISMATCH,
                         feature_group = ndfdata$MATCH_INDEX,
                         feature_class = ndfdata$PROBE_CLASS)

  ndfdata <- ndfdata2
  rm(ndfdata2)

  # Getting only the upper left corner
  ndfdata <- ndfdata[order(ndfdata$feature_id,ndfdata$Y,ndfdata$X),]
  dup <- duplicated(ndfdata$feature_id)
  ndfdata <- ndfdata[!dup,]

  # Putting in the same order as the XYS
  ndfdata <- ndfdata[order(ndfdata$Y,ndfdata$X),]
  return(ndfdata)
}
