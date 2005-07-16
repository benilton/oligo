## BC: Jul 13, "pd" at the begining
cleanPlatformName <- function(x)
  gsub("[_-]","",paste("pd",tolower(x),sep=""))

## BC: added nrow/ncol
setClass("platformDesign",
         representation(featureInfo = "environment",
                        manufacturer = "character",
                        type = "character",
                        nrow = "numeric",
                        ncol = "numeric"))

##functions and methods
featureInfo <- function(object) object@featureInfo

## i just made names a generic method... hope this isnt bad!
if( is.null(getGeneric("names")))
  setGeneric("names", function(x)
             standardGeneric("names"))
setMethod("names","platformDesign",
          function(x){
            return(ls(featureInfo(x)))
          })

setMethod("[", "platformDesign", function(x, i, j, ..., drop=FALSE) {
  if( missing(j) ) j <- names(x)  else j <- names(x)[j] 
  if( missing(i) ) i <- 1:nProbes(x)

  env <- new.env()

  for(col in j){
    val <- get(col,featureInfo(x))[i]
    assign(col,val,col,envir=env)
  }
  
  return(new("platformDesign",featureInfo=env,manufacturer=x@manufacturer,type=x@type))
})

##show method
setMethod("show","platformDesign", function(object){
  cat("Manufacturer:",object@manufacturer,"\n")
  cat("Array type:",object@type,"\n")
  cat("Number of features:",nProbes(object),"\n")
})


as.data.frame.platformDesign <- function(x,row.names=NULL,optional=FALSE){

  myCall <- "df=data.frame("
  for(i in seq(along=names(x))){
    myCall <- paste(myCall,"get(\"",names(x)[i],"\",featureInfo(x))",sep="")
    if(i==length(names(x))){myCall <- paste(myCall,")")} else{ myCall<-paste(myCall,",")}
  }
  eval(parse(text=myCall))
  names(df) <- names(x)
  df <- as.data.frame(df,row.names=row.names,optional=optional)
  return(df)
}

"$.platformDesign" <- function(object, val)
  get(val,featureInfo(object))
      


##X is arbitrarily chosen
nProbes <- function(object){
  colname <- ls(featureInfo(object))[1]
  length( get(colname,featureInfo(object)))
}

##probeNames method
if( is.null(getGeneric("probeNames")))
  setGeneric("probeNames", function(object, ...)
             standardGeneric("probeNames"))

setMethod("probeNames","platformDesign",
          function(object){
            return(get("feature_name",featureInfo(object)))
          })
                                   

##pmindex method
if( is.null(getGeneric("pmindex")))
  setGeneric("pmindex", function(object,...)
             standardGeneric("pmindex"))

###NOTE: THIS WILL CHANGE CAUSE feature type will be a vector
setMethod("pmindex", "platformDesign",
          function(object){
            Index=which(get("feature_type",featureInfo(object))=="PM")
            pns = probeNames(object)
            return(Index[order(pns[Index])])
          })

