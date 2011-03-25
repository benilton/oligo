###########################################################
##
## file: fitPLM.R
##
## Copyright (C) 2003-2008   Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Jan 15, 2003
##
## Last modified: See History below
##
## method for fitting robust linear models.
## interface to C code which does the heavy lifting
##
##
## by default model should be written in the following terms
##
## PM  : should appear on left hand side
## ~ : equivalent to  =
## -1 : remove the intercept term
## probes : include terms for probes in the model. a variable of type factor
## samples : defaults to be a factor, one for each chip or array
## other names for chip level factors or covariates. By default we will check the
## phenoData object to see if the named covariate or factor exists in the
## phenoData slot of the provided affybatch, then the parent enivronment is checked
## if the named variable does not exist there then procedure exits.
##
## an example model fitting terms for probes and samples
##
##    pm ~ -1 + probes + samples
##
## another example model
##
##    pm ~ -1 + probes + trt.group
##
## where trt.group is a factor in the phenoData object.
##
## A valid model should always include some chip effects
##
##
## Feb 1, 2003 - cache rownames(PM(object)) since this is slow operation
## Feb 15 - set model description to the function call: THIS was undone.
##         made it possible to fit models with an intercept term.
## Feb 17 -  start testing and fixing up column naming when model fit with no
##         intercept term, also make sure that se.exprs is set if there is
##         an intercept in the fitted model. re-order some parts of the function
##         work towards being able to fit models with no slope or probe-effects.
## Feb 22 - Add ability to fit models where parameters use a sum to zero constraint.
##          A similar mechanism to the variables.type parameter will be used.
##          note that R uses the convention for endpoint constraint that first
##          item is constrained to 0 and contr.sum that the last item is left unshown
##          we will try to follow this convention, except in the case constraints on
##          the probe coefficient.
##          Remember that the probe coefficient is handled by c code.
## Mar 22 - Add in ability to use LESN background correction.
## Mar 23 - Add ability to use MAS style background in fit PLM
## Jun 4-8 - Ability to specify other psi functions.
## Jul 22 - No longer need to specify psi.k, if left unspecified (ie null)
##          we set it to the default parameter. These defaults are the following:
## ** Huber - k = 1.345
## ** Fair - k = 1.3998
## ** Cauchy - k=2.3849 
## ** Welsch - k = 2.9846
## ** Tukey Biweight - k = 4.6851
## ** Andrews Sine - K = 1.339
##          the other methods have no tuning parameter so set to 1.
## Jul 27 - Clean up parameters passed to function
## Sep 2  - Residuals are now stored in the outputed PLMset
## Sep 4  - may toggle what is actually outputed
##          in particular in respect to weights, residuals
##          var.cov and resid.SE
## Sep6 - Sept 8  -  More changes to how naming is done.
## Sep 12 - Fix how constraints are applied when factor variable
##          is defined in parent enivronment. Basically the constraint
##          was being ignored. Note that this fix was also applied to
##          the 0.5-14 Series package.
## Sep 12 - model.param was introduced as argument
##          se.type, psi.type, psi.k were placed into
##          this parameter
## Oct 10 - fix spelling of McClure
## Jan 18, 2004 - move some of the sub functions to internalfunctions.R
## Feb 23, 2004 - subset paramter added
## Mar 1,  2004 - moved chip-level design matrix creation to its own function
## Mar 2, 2004 - rewrote design matrix function to handle interaction terms 
## May 26, 2004 - allow use of MM as response term in fitPLM
## May 27, 2004 - if default model ie -1 + probes + samples flag that a optimized
##                algorithm should be used.
## Jun 28, 2004 - allow MM or PM as a covariate term in the model
## Jul 14, 2004 - introduce PLM.designmatrix3 for dealing with new fitting algorithms
## Aug 3, 2004 -  a new fitPLM introduced along with various support functions.
## Feb 18, 2005 - gcrma background is also now an option
## Apr 27-28, 2006 - clean up how normalization methods are check and normalization parameters are validated.
## Oct 10, 2006 - add verbosity.level argument to fitPLM
## Jan 3, 2007 - exprs, and se.exprs are no longer used in PLMset
## Aug 10, 2007 - replace error with stop
## Dec 1, 2007 - comment out fitPLM.old, PLM.designmatrix, PLM.designmatrix2 (will be removed in later version and all are long defunct)
## Feb 3, 2007 - remove all commented out code. Make sure fitPLM stored narrays in PLMset object
## Jan 20, 2009 - Fix issue where factor variable not coerced to integer
## 
###########################################################





###
###
### check the model and supplied constraint.type vector (or list)
### output a constraint.type vector 
###

verify.constraint.types <- function(model,constraint.type){

  if (!(is.vector(constraint.type) | is.list(constraint.type))){
    stop("constraint.type must be a vector or list.")
  }
  model.terms <- terms(model)
  ct.type <- NULL
  if (!any(is.element(c("Default","default"),names(constraint.type)))){
    cat("Warning: No default constraint specified. Assuming 'contr.treatment'.\n")
    ct.type <- c(ct.type,"default"="contr.treatment")
  }

  if (length(names(constraint.type)) != length(constraint.type)){
    stop("constraint.type is incorrectly named\n")
  }

  if (any(is.element("",names(constraint.type)))){
    stop("Some elements are not named in constraint.type")
  }

  if (any(!is.element(constraint.type,c("contr.treatment","contr.sum")))){
    stop("Constraints must be either 'contr.treatment' or 'contr.sum'")
  }
  

  if (is.list(constraint.type)){
    constraint.type <- do.call(c,constraint.type)
  }
  

  return(c(ct.type,constraint.type))

}



##  verify.constraint.types(PM ~ -1 + probes + samples,c(blah="blah"))


verify.variable.types <- function(model,variable.type){

  vt <- NULL

  if (!(is.vector(variable.type) | is.list(variable.type))){
    stop("variable.type must be a vector or list.")
  }
  
  if (!any(is.element(c("Default","default"),names(variable.type)))){
    cat("Warning: No default variable type so assuming 'factor'\n")
    vt <- c(vt,"default"="factor")
  }

  if (length(names(variable.type)) != length(variable.type)){
    stop("variable.type is incorrectly named\n")
  }


  if (is.list(variable.type)){
    variable.type <- do.call(c,variable.type)
  }

  if (any(is.element("",names(variable.type)))){
    stop("Some elements are not named in variable.type")
  }
  
  if (any(!is.element(variable.type,c("factor","covariate")))){
    stop("variable.type must be either 'factor' or 'covariate'")
  }

  for (cur.vt in c("probes","probe.type","samples")){
    if (is.element(cur.vt,names(variable.type))){
      if (variable.type[cur.vt] == "covariate"){
        stop(cur.vt, " must not be specified as being a 'covariate'. It must be treated as 'factor'.")
      }
    }
  }

  for (cur.vt in c("mm","MM","pm","PM"))
    if (is.element(cur.vt,names(variable.type))){
      if (variable.type[cur.vt] == "factor")
        stop(cur.vt," variable must be a 'covariate'.")
    }

  
  return(c(vt,variable.type))

}



#verify.variable.types(PM ~ -1 + probes + samples,c(blah="blah"))



verify.model.param <- function(object,model,subset=NULL,model.param=list()){

  defaults <- list(trans.fn="log2", se.type = 4, psi.type = 0, psi.k =NULL,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL)

  if (!is.list(model.param)){
    stop("model.param must be a list")
  }

  if (length(model.param)==0){
    if (is.null(defaults["psi.k"][[1]])){
      defaults["psi.k"] <- get.default.psi.k(defaults["psi.type"][[1]])
    }
    return (defaults)
  }

  if (any(!is.element(names(model.param),c("trans.fn","se.type","psi.type","psi.k","max.its","init.method","weights.chip","weights.probe")))){
    stop("model.param should only have items named: ,",paste("trans.fn","se.type","psi.type","psi.k","max.its","init.method","weights.chip","weights.probe",sep=", "))
  }

  
  for (item in names(model.param)){
    if (item == "trans.fn"){
      if (!is.element(model.param[item][[1]],c("log2","loge","ln","log10","sqrt","cuberoot"))){
        stop("trans.fn in model.param should be one of: ",paste("log2","loge","ln","log10","sqrt","cuberoot",sep=", "))
      } else {
        defaults[item] <- model.param[item]
      }
    }
    if (item == "se.type"){
      if (!is.numeric(model.param[item][[1]])){
        stop("se.type should be numeric")
      } else if (!is.element(model.param[item][[1]],1:4)){
        stop("se.type should be 1,2,3,4")
      } else {
        defaults[item] <- model.param[item]
      }
    }
    if (item == "psi.type"){
      if (is.numeric(model.param[item][[1]])){
        if (!is.element(model.param[item][[1]],0:6)){
          stop("psi.type should be a string")
        } else {
          defaults[item] <- model.param[item]
        }
      } else {
        defaults[item] <- get.psi.code(model.param[item])
      }

    }

    if (item == "max.its"){
      if (!is.numeric(model.param[item][[1]])){
        stop(item," should be numeric")
      } else {
        defaults[item] <- model.param[item]
      }
    }
    if (item == "init.method"){
      if (!is.element(model.param[item][[1]],c("ls","Huber"))){
        stop("init.method should be one of: ls, Huber")
      } else {
        defaults[item] <- model.param[item]
      }
    }
    if (item == "weights.chip"){
      if (!is.null(model.param[item][[1]])){
        if (!is.vector(model.param[item][[1]])){
          stop("weights.chip should be a vector")
        } else if (length(model.param[item][[1]]) != dim(exprs(object))[2]){
          stop("weights.chip is of incorrect length")
        } else {
          defaults[item] <- model.param[item]
        }
      }
    }
    if (item == "weights.probe"){
      if (!is.null(model.param[item][[1]])){
        response.term <- attr(terms(model),"variables")[[2]]
        if (is.element(as.character(response.term),c("MM","PM","pm","mm"))){
          weights.probe.length <- dim(pm(object,subset))[1]
        } else {
          weights.probe.length <- 2*dim(pm(object,subset))[1]
        }
        if (!is.vector(model.param[item][[1]])){
          stop("weights.probe should be a vector")
        } else if (length(model.param[item][[1]]) !=weights.probe.length) {
          stop("weights.probe is of incorrect length")
        } else {
           defaults[item] <- model.param[item]
        }
      }
    }
  }
  
  if (!is.element("psi.k",names(model.param))){
    if (is.null(defaults["psi.k"][[1]])){
      defaults["psi.k"] <- get.default.psi.k(defaults["psi.type"][[1]])
    }
  } else {
    defaults["psi.k"] <- model.param["psi.k"]
  }
  
  return(defaults)


}

##verify.model.param(list(trans.fn="cuberoot", se.type = 4, psi.type = 0, psi.k =1.345,max.its = 20, init.method = "ls",weights.chip=NULL,weights.probe=NULL))








verify.output.param <- function(output.param=list()){
  
  if (!is.list(output.param)){
    stop("output.param must be a list")
  }

  if (length(output.param) == 0){
    ### this is the default output
    return(list(weights = TRUE, residuals = TRUE, varcov ="none", resid.SE = TRUE))
  }

  if (any(!is.element(names(output.param),c("weights","residuals","varcov","resid.SE")))){
    stop("Only items named 'weights', 'residuals', 'varcov' or 'resid.SE' should be in output.param")
  }

  out.list <- list(weights = TRUE, residuals = TRUE, varcov ="none", resid.SE = TRUE)

  for (item in names(output.param)){
    if (item == "weights"){
      if (is.logical(output.param[item][[1]])){
        out.list[item] <- output.param[item]
      } else {
        stop("Logical needed for ",item," in output.param.")
      }
    }
    if (item == "residuals"){
      if (is.logical(output.param[item][[1]])){
        out.list[item] <- output.param[item]
      } else {
        stop("Logical needed for ",item," in output.param.")
      }
    }
    if (item == "varcov"){
      if (!is.element(output.param[item][[1]],c("none","chiplevel","all"))){
        stop("varcov in ouput.param must be: 'none', 'chiplevel' or 'all'")
      } else {
        out.list[item] <- output.param[item]
      }
    }
    if (item == "resid.SE"){
      if (is.logical(output.param[item][[1]])){
        out.list[item] <- output.param[item]
      } else {
        stop("Logical needed for ",item," in output.param.")
      }
    }   
  }
  return(out.list)

}


###
### PLM.designmatrix3  -  Third generation function. Produces  a list
### 
###  
###
### Input: the AffyBatch
###        a model formula
###        a variable types list
###        a constraints list
###
### Output: a list which must contain all of the following items
###             mmorpm.covariate  :  an integer  either -1  meaning PM is covariate, 0 means no covariate,   1 meaning MM is covariate
###             response.variable :  an integer  either -1  meaning MM is response,  0 means PM and MM are response,  1 is PM response variable
###             which.parameter.types : an integer vector of length 5 values should be only 0 or 1
###                                    element 1 - 0 means no intercept, 1 means intercept
###                                    element 2 - 0 means no chip-level factor or covariate variables, 1 means there are. Only one of elements 2 and 3 should have value 1
###                                    element 3 - 0 means no sample effects, 1 means there are sample effects
###                                    element 4 - 0 means no probe.type effect,  1 means probe.type effect
###                                    element 5 - 0 means no probe effect, 1 means probe effect
###             strata: an integer vector of length 5,
###                                    element 1 - ignored
###                                    element 2 - ignored
###                                    element 3 - ignored
###                                    element 4 - 0 means overall probe.type effect,
###                                                1 means sample specific probe.type effect,
###                                                2 means within levels of a treatment/genotype factor variable
###                                                OTHER INVALID
###                                    element 5 - 0 means overall probe effects
###                                                1 INVALID
###                                                2 means within levels of a treatment/genotype factor variable
###                                                3 means within probe.types
###                                                4 means within probe.types within levels of treatment/genotype factor variable
###                                                OTHER INVALID
###              constraints: an integer vector of length 5,
###                                    element 1 - ignored
###                                    element 2 - ignored
###                                    element 3 - sample effect: 0 means unconstrainted, -1 means sum to zero, 1 means contr.treatment
###                                                OTHER INVALID
###                                    element 4 - probe.type effect: 0 means unconstrainted, -1 means sum to zero, 1 means contr.treatment
###                                                OTHER INVALID
###                                    element 5 - probe effect:  0 means unconstrainted, -1 means sum to zero, 1 means contr.treatment
###                                                OTHER INVALID
###              probe.type.trt.factor:  a vector of the same length as the number of arrays.
###                                      values should be between 0 and max.probe.type.trt.factor
###                                      the values should correspond to which level of the treatment or genotype variable is assigned to the respective array     
###              max.probe.type.trt.factor: the maximum value in "probe.type.trt.factor"
###              probe.type.levels: A list of only one element. The element should be named for the factor and it should contain a vector of character strings
###                                of length max.probe.type.trt.factor with the names being those of the levels of the factor
###              probe.trt.factor: a vector of the same length as the number of arrays.
###                                values should be between 0 and max.probe.trt.factor
###                                the values should correspond to which level of the treatment or genotype variable is assigned to the respective array
###              max.probe.trt.factor: the maximum value in "probe.trt.factor"
###              probe.trt.levels: A list of only one element. The element should be named for the factor and it should contain a vector of character strings
###                                of length max.probe.trt.factor with the names being those of the levels of the factor
###              chipcovariates: a matrix.
###                              if which.parameter.types[2] == 0 then should be matrix(0,0,0)
###                              otherwise it should be a matrix containing values of chip-level factor variables
###                              for each array (ie dim n_arrays * n_chiplevelcovariates)
###
### There are other rules governing what is and what is not a permissible model. Please see other documentation for those rules.

PLM.designmatrix3 <- function(object,model=PM ~ -1 + probes + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment")){

  ## a function for working out which constraints are applied to a named term
  ## note that "probes", "probe.type" have
  ## special defaults that will used irrespective of
  ## default unless there is a specifically specified constraint
  ## for that parameter.
  ## probes -- "contr.sum"
  ## probe.type -- "contr.treatment"
  ##
  ##
  which.constraint <- function(term,constraint.type){
    if (term == "probes"){
      if (is.element("probes",names(constraint.type))){
        if (constraint.type[term]=="contr.sum"){
          return("contr.sum")
        } else if (constraint.type[term] =="contr.treatment"){
          return("contr.treatment")
        } else {
          stop("The constraint type ",names(constraint.type)[term]," for ", term, "is not understood.")
        }
      } else {
        return("contr.sum")
      }
      
    }
    if (term == "probe.type"){
      if (is.element("probe.type",names(constraint.type))){
        if (constraint.type[term]=="contr.sum"){
          return("contr.sum")
        } else if (constraint.type[term] =="contr.treatment"){
          return("contr.treatment")
        } else {
          stop("The constraint type ",names(constraint.type)[term]," for ", term, "is not understood.")
        }
      } else {
        return("contr.treatment")
      }
    }
    
    if (is.element(term,names(constraint.type))){
      if (constraint.type[term]=="contr.sum"){
        return("contr.sum")
      } else if (constraint.type[term] =="contr.treatment"){
        return("contr.treatment")
      } else {
        stop("The constraint type ",names(constraint.type)[term]," for ", term, "is not understood.")
      }
    } else {
      if (any(is.element(names(constraint.type),"default"))){
        return(constraint.type["default"])
      } else {
        warning("No default constraint was supplied so assuming : contr.treatment for ",term)
        return("contr.treatment")
      }
    }
  }
  which.variabletype <- function(term,variable.type){
    if (is.element(term,names(variable.type))){
      if (variable.type[term]=="covariate"){
        return("covariate")
      } else if (variable.type[term] =="factor"){
        return("factor")
      } else {
        stop("The variable type ",names(variable.type)[term]," for ", term, "is not understood.")
      }
    } else {
      if (any(is.element(names(variable.type),"default"))){
        return(variable.type["default"])
      } else {
        warning("No default variable type was supplied so assuming : factor for ",term)
        return("factor")
      }
    }
  }
  


  
  model.terms <- terms(model)
  mt.variables <- attr(model.terms,"variables")

  ## establish the response variable
  if (!is.element(as.character(mt.variables[[2]]),c("PM","pm","MM","mm","PMMM","pmmm"))){
    stop(paste("Response term in model should be 'PM' or 'pm' or 'MM' or 'mm' or 'PMMM' or 'pmmm'."))
  }
  
  if (is.element(as.character(mt.variables[[2]]),c("PM","pm"))){
    response.variable <- 1
  } else if (is.element(as.character(mt.variables[[2]]),c("MM","mm"))){
    response.variable <- -1  
  } else   if (is.element(as.character(mt.variables[[2]]),c("PMMM","pmmm"))){
    response.variable <- 0
  } else {
    stop("Response term in model should be 'PM' or 'pm' or 'MM' or 'mm' or 'PMMM' or 'pmmm'.")
  }

  ## check to see if there is a PM or MM covariate

  mt.termlabels <- attr(model.terms,"term.labels")

  if (any(is.element(as.character(mt.termlabels),c("PM","pm")))){
    mmorpm.covariate <- -1
  } else if (any(is.element(as.character(mt.termlabels),c("MM","mm")))){
    mmorpm.covariate <- 1
  } else {
    mmorpm.covariate <- 0
  }

  ## check to see if there is a PMMM covariate (which is not allowable)
  if (any(is.element(as.character(mt.termlabels),c("PMMM","pmmm")))){
     stop("Cannot have PMMM or pmmm as a covariate variable.")
   }

  
  ## check to see if there is a PM or MM covariate when response is PMMM (note we will allow PM ~ PM and MM ~ MM stupidity)

  if (response.variable ==0 & mmorpm.covariate!=0){
    stop("Cannot have PMMM as response and MM or PM as a covariate.")
  }



  

  ## initialize some of the output
  which.parameter.types <- rep(0,5)
  strata <- rep(0,5)
  constraints <- rep(0,5)

  probe.type.trt.factor <- rep(0,dim(intensity(object))[2])
##  probe.type.trt.factor <- rep(0,dim(exprs(object))[2])
  max.probe.type.trt.factor <- 0
  probe.type.levels <- list()
  probe.trt.factor <- rep(0,dim(intensity(object))[2])
##  probe.trt.factor <- rep(0,dim(exprs(object))[2])
  max.probe.trt.factor <- 0
  probe.trt.levels <- list()
  ## now check to see what else we have in the model
  ## intercept??
  
  mt.intercept <- attr(model.terms,"intercept")
  if (mt.intercept){
    which.parameter.types[1] <- 1
  }
  
  mt.variables <- as.list(mt.variables)[3:length(mt.variables)]  

  ## check to see if probe.type is in the model when PMMM is not response

  if (response.variable != 0 & any(is.element(as.character(mt.variables),"probe.type"))){
    stop("Cannot have 'probe.type' without PMMM response variable.")
  }

  has.nonspecialvariables <- TRUE

  if (has.nonspecialvariables){
    ## check to see if the nonspecialvariables exist
    for (variable in mt.variables){
      variable <- as.character(variable)
      if (is.element(variable,c("PM","pm","MM","mm","probes","probe.type","samples"))){
        next
      }  
      ## check if variable is in phenoData
      if (is.element(variable,  names(pData(object)))){
        next;
      }
      ## check to see if variable exists elsewhere
      if (exists(variable)){
        # check it is of appropriate length
        if (length(eval(as.name(variable))) !=dim(intensity(object))[2]){
##        if (length(eval(as.name(variable))) !=dim(exprs(object))[2]){
          stop("The variable ",variable," is of the wrong length.")
        }
        next;
      }
      stop("The variable ",variable," does not seem to exist.")
    }

    
    mt.factors <- attr(model.terms,"factors")
    mt.order <- attr(model.terms,"order")
    if (mt.intercept){
      the.formula <- "~"
      firstfactor <- FALSE
    } else {
      the.formula <- "~ -1 +"
      firstfactor <- TRUE
    }
    the.frame <- NULL
    the.frame.names <- NULL
    terms.names <- NULL
    ## lets make the model matrix for chipcovariates and any other interaction with probes, probe.types
    i <- 0
    for (term in colnames(mt.factors)){
      i <- i +1
     # if (is.element(term ,c("PM","pm","MM","mm","probes","probe.type","samples"))){
        ## special variable so do nothing
    ##    next
     # }
      if (mt.order[i] ==1){

        if (is.element(term,c("PM","pm","MM","mm","probes","probe.type","samples"))){
          ## special variable so do nothing but set the correct flag
          if (is.element(term,"samples")){
            if (which.parameter.types[3]){
              stop("Can't have 'samples' appear in more than one term in model.")
            }
            which.parameter.types[3] <- 1
            ct <- which.constraint("samples",constraint.type)
            if (ct == "contr.sum"){
              constraints[3] <- -1
            } else {
              constraints[3] <- 1
            }
          }
          if (is.element(term,"probe.type")){
            if (which.parameter.types[4]){
               stop("Can't have 'probe.type' appear in more than one term in model, except with probes")
             }
            which.parameter.types[4] <- 1
            ct <- which.constraint("probe.type",constraint.type)
            if (ct == "contr.sum"){
              constraints[4] <- -1
            } else {
              constraints[4] <- 1
            }

            
          }
          if (is.element(term,"probes")){
            if (which.parameter.types[5]){
              stop("Can't have 'probes' appear in more than one term in model.")
            }
            which.parameter.types[5] <- 1
            ct <- which.constraint("probes",constraint.type)
            if (ct == "contr.sum"){
              constraints[5] <- -1
            } else {
              constraints[5] <- 1
            }
          }
          next
        }
        if (is.element(term,  names(pData(object)))){
          the.frame <- cbind(the.frame,pData(object)[,names(pData(object))==term])
          trt.values <- pData(object)[,names(pData(object))==term]
        } else {
          the.frame <- cbind(the.frame,eval(as.name(term)))
          trt.values <- eval(as.name(term))
        }
        
        
        if (which.variabletype(term,variable.type)=="covariate"){
          # this variable is a covariate
          the.formula <- paste(the.formula,term)
          terms.names <- c(terms.names,term)
        } else {
          ##this variable is a factor variable
          this.constraint <- which.constraint(term,constraint.type)
          if (this.constraint == "contr.treatment"){
            the.formula <- paste(the.formula,"+","C(as.factor(",term,"),","contr.treatment",")")
            if (!firstfactor){
              for (levs in levels(as.factor(trt.values))[-1]){
                terms.names <- c(terms.names,paste(term,"_",levs,sep=""))
              }
            } else {
              for (levs in levels(as.factor(trt.values))){
                terms.names <- c(terms.names,paste(term,"_",levs,sep=""))
              }
              firstfactor <- FALSE
            }           
          } else {
            the.formula <- paste(the.formula,"+", "C(as.factor(",term,"),","contr.sum",")")
            if (!firstfactor){
              for (levs in levels(as.factor(trt.values))[-length(levels(as.factor(trt.values)))]){
                terms.names <- c(terms.names,paste(term,"_",levs,sep=""))
              }
            } else {
              for (levs in levels(as.factor(trt.values))){
                terms.names <- c(terms.names,paste(term,"_",levs,sep=""))
              }
              firstfactor <- FALSE
            }
          }
        }
        the.frame.names <- c(the.frame.names,term)
      } else {
        # a higher term 
        in.term <- strsplit(term,":")[[1]]

        if (any(is.element(in.term,c("probes","probe.type")))){
          ## a so called special case need to handle appropriately
          if (any(is.element(in.term,"probes"))){
            # figure out if probes is within a treatment variable, probe.type or both.
            if (which.parameter.types[5]){
              stop("Can't have 'probes' appear in more than one term in model.")
            }
            which.parameter.types[5] <- 1
            in.term <- in.term[!is.element(in.term,"probes")]
            if (length(in.term) > 2){
              stop("Can't estimate probe effect within more than two variables.")
            }
            if (length(in.term) ==2){
              if (!any(is.element(in.term,"probe.type"))){
                stop("Can't estimate probe effect within two variables without one being probe.type,") 
              }
              in.term <- in.term[!is.element(in.term,"probe.type")]
              if (which.variabletype(in.term,variable.type) == "covariate"){
                stop("Can't have interaction terms involving covariate variables.")
              }
              
              if (is.element(in.term,  names(pData(object)))){
                in.levels <- as.factor(pData(object)[,names(pData(object)) == in.term])
              } else {
                in.levels <- as.factor(eval(as.name(in.term)))
              }
              ## check that there is at least two arrays
              if (any(tabulate(in.levels) < 2)){
                stop("Need to have at least two arrays for each level of ",in.term, " when estimating probes within levels of this variable." )
              }
              
              # now build what we need for the output
              max.probe.trt.factor <- nlevels(in.levels)-1
              probe.trt.factor <- NULL
	      in.levels.2 <- as.numeric(in.levels)
              for (level in in.levels.2){
                probe.trt.factor <- c(probe.trt.factor,level-1)
              }
              probe.trt.levels <- list(in.term=levels(in.levels))
              names(probe.trt.levels) <- in.term
              strata[5] <- 4
              ct <- which.constraint("probes",constraint.type)
              if (ct == "contr.sum"){
                constraints[5] <- -1
              } else {
                constraints[5] <- 1
              }
            } else if (length(in.term) ==1){
              if (is.element(in.term,"probe.type")){
                strata[5] <- 3
                ct <- which.constraint("probes",constraint.type)
                if (ct == "contr.sum"){
                  constraints[5] <- -1
                } else {
                  constraints[5] <- 1
                }
              } else {
                if (which.variabletype(in.term,variable.type) == "covariate"){
                  stop("Can't have interaction terms involving covariate variables.")
                }
                
                if (is.element(in.term,  names(pData(object)))){
                  in.levels <- as.factor(pData(object)[,names(pData(object)) == in.term])
                } else {
                  in.levels <- as.factor(eval(as.name(in.term)))
                }

                if (any(tabulate(in.levels) < 2)){
                  stop("Need to have at least two arrays for each level of ",in.term, " when estimating probes within levels of this variable." )
                }
                max.probe.trt.factor <- nlevels(in.levels)-1
                probe.trt.factor <- NULL
		in.levels.2 <- as.numeric(in.levels)
                for (level in in.levels.2){
                  probe.trt.factor <- c(probe.trt.factor,level-1)
                }
                probe.trt.levels <- list(levels(in.levels))
                names(probe.trt.levels) <- in.term
                strata[5] <- 2
                ct <- which.constraint("probes",constraint.type)
                if (ct == "contr.sum"){
                  constraints[5] <- -1
                } else {
                  constraints[5] <- 1
                }
              }
            }
          } else if (any(is.element(in.term,"probe.type"))){
            if (which.parameter.types[4]){
              stop("Can't have 'probe.type' appear in more than one term in model, except with probes")
            }
            which.parameter.types[4] <- 1
            ## figure out if probe.type is in samples or a treatment factor
            in.term <- in.term[!is.element(in.term,"probe.type")]
            if (length(in.term) > 1){
              stop("Can't estimate probe.type within more than one variable.")
            }
            if (is.element(in.term,"samples")){
              strata[4] <- 1
              ct <- which.constraint("probe.type",constraint.type)
              if (ct == "contr.sum"){
                constraints[4] <- -1
              } else {
                constraints[4] <- 1
              }
            } else{
              strata[4] <- 2
              ct <- which.constraint("probe.type",constraint.type)
              if (ct == "contr.sum"){
                constraints[4] <- -1
              } else {
                constraints[4] <- 1
              }
              if (which.variabletype(in.term,variable.type) == "covariate"){
                stop("Can't have interaction terms involving covariate variables.")
              }
              if (is.element(in.term,  names(pData(object)))){
                in.levels <- as.factor(pData(object)[,names(pData(object)) == in.term])
              } else {
                in.levels <- as.factor(eval(as.name(in.term)))
              }
              max.probe.type.trt.factor <- nlevels(in.levels)-1
              probe.type.levels <- list(levels(in.levels))
              names(probe.type.levels) <- in.term
              probe.type.trt.factor <- NULL
	      in.levels.2 <- as.numeric(in.levels)
              for (level in in.levels.2){
                probe.type.trt.factor <- c(probe.type.trt.factor,level-1)
              }
            }
          } 
        } else {
          ## a more general interaction term
          if (any(is.element(in.term,c("samples","PM","MM","pm","mm")))){
            stop("Cannot have 'samples','PM' or 'MM' in an interaction term.")
          }
          for (current.term in in.term){
            if (which.variabletype(term,variable.type)=="covariate")
              stop("Can't have interaction terms involving covariate variables.")
          }
          for (current.term in in.term){
            if (is.element(current.term,  names(pData(object)))){
              the.frame <- cbind(the.frame,pData(object)[,names(pData(object))==current.term])
            } else {
              the.frame <- cbind(the.frame,eval(as.name(current.term)))
            }
            the.frame.names <- c(the.frame.names,current.term)
          }
          
          ## now lets make the actual formula   
          this.term <- "C("
          for (current.term in in.term){
            this.term <- paste(this.term,"as.factor(",current.term,")")
            if (current.term != in.term[length(in.term)]){
              this.term <- paste(this.term,":")
            }
          }

          if (which.constraint(term,constraint.type) == "contr.treatment"){
            this.term <- paste(this.term,",contr.treatment)")
          } else {
            this.term <- paste(this.term,",contr.sum)")
          }
          #terms.names <- c(terms.names,interaction.terms)
          


          firstterm <- TRUE
          interaction.terms <- as.vector("")
          for (current.term in in.term){
            if (is.element(current.term,  names(pData(object)))){
              trt.values <- pData(object)[,names(pData(object))==current.term]
            } else {
              trt.values <- eval(as.name(current.term))
            }
            
            levs <- levels(as.factor(trt.values))
            if (!firstterm){
              interaction.terms <-  as.vector(t(outer(interaction.terms,paste(":",current.term,"_",levs,sep=""),paste,sep="")))
            } else {
              interaction.terms <-  as.vector(t(outer(interaction.terms,paste(current.term,"_",levs,sep=""),paste,sep="")))
              firstterm <- FALSE
            }
          }
          if (!firstfactor){
            if (which.constraint(term,constraint.type) == "contr.treatment"){
              interaction.terms <- interaction.terms[-1]
            } else {
              interaction.terms <- interaction.terms[-length(interaction.terms)]
            }
          } else {
            firstfactor <- FALSE
          }
          
          terms.names <- c(terms.names,interaction.terms)    
                    
          the.formula <- paste(the.formula,"+",this.term)
        }
      }
    }

    if (!is.null(the.frame)){ 
      the.frame <- as.data.frame(the.frame)
      colnames(the.frame) <- the.frame.names  
      
      ##print(the.formula)
      ##print(the.frame)
                                        #print(as.list(ct))
                                        #print(model.matrix(as.formula(the.formula),the.frame))
      chip.covariates <- model.matrix(as.formula(the.formula),the.frame)

        
      if (qr(chip.covariates)$rank < ncol(chip.covariates)){
        stop("chip-level effects appear to be singular: singular fits are not implemented in fitPLM")
      }
  
      if (mt.intercept){
        chip.covariates <- as.matrix(chip.covariates[,-1])
      }
      colnames(chip.covariates) <- terms.names
      
      which.parameter.types[2] <- 1
    } else {
      chip.covariates <- matrix(0,0,0)
    }
  } else {
    chip.covariates <- matrix(0,0,0)
  }


  if ((which.parameter.types[2] == 1) & (which.parameter.types[3] == 1)){
    stop("Can't have 'samples' and chip-level variables in same model.")
  }
  
  # final checks for constraints, unconstrain using the order intercept:chip-level:samples:probe.types:probes

  if (all(which.parameter.types[1:2] == 0)){
    constraints[3] <- 0
  }

  if (all(which.parameter.types[1:3] == 0)){
    constraints[4] <- 0
  }

  if (all(which.parameter.types[1:4] == 0)){
    constraints[5] <- 0
  }

  
  list(mmorpm.covariate=mmorpm.covariate,response.variable=response.variable,which.parameter.types=as.integer(which.parameter.types),strata=as.integer(strata),constraints=as.integer(constraints),probe.type.trt.factor=as.integer(probe.type.trt.factor),max.probe.type.trt.factor=max.probe.type.trt.factor,probe.type.levels=probe.type.levels,probe.trt.factor=as.integer(probe.trt.factor),max.probe.trt.factor=max.probe.trt.factor,probe.trt.levels=probe.trt.levels,chipcovariates =chip.covariates)
  


}


#PLM.designmatrix3(Dilution,model=MM ~ PM -1 + probes+ liver+scanner)
#PLM.designmatrix3(Dilution,model=MM ~ PM -1 + probes+ probe.type:liver)

verify.bg.param <- function(R.model, background.method,background.param = list()){

  bg.code <- get.background.code(background.method)
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}
  LESN.param <-list(baseline=0.25,theta=4)
  LESN.param <- convert.LESN.param(LESN.param)
  
  default.b.param <- list(type="separate",densfun =  body(bg.dens), rho = new.env(),lesnparam=LESN.param,ideal=NULL)

  if (R.model$response.variable == 0){
    response.variable <- "PMMM"
  } else if (R.model$response.variable == 1){
    response.variable <- "PM"
    if (R.model$mmorpm.covariate == 0){
      default.b.param["type"] <- "pmonly"
    } 
  } else if (R.model$response.variable == -1){
    response.variable <- "MM"
    if (R.model$mmorpm.covariate == 0){
      default.b.param["type"] <- "mmonly"
    }
  }

  if (is.element(as.character(response.variable),c("PMMM","pmmm"))){
    if (length(names(background.param)) !=0){
      if (is.element(names(background.param),"type")){
        if (is.element(background.param["type"],c("pmonly","mmonly"))){
          stop("You can not apply a background 'pmonly' or 'mmonly' with the 'PMMM' response.")
        } else {
          if (is.element(background.param["type"],c("separate","together"))){
            default.b.param["type"] <- background.param["type"]
          } else {
            stop("type should be 'separate','together','pmonly' or 'mmonly'.")
          }
        }
      }
    }
    if (is.element(background.method,c("IdealMM","MASIM"))){
      stop("Can't use the Ideal Mismatch background corrections with the 'PMMM' response.") 
    }
  }

  if (!all(is.element(names(background.param),c("type","lesnparam")))){
    stop("Unknown parameters in background.param")
  }

  
  if (is.element(background.param["type"],c("separate","pmonly","mmonly")) & is.element(background.method,c("MAS","MASIM","IdealMM"))){
    warning("'together' type has been used in place of 'separate', 'pmonly' or 'mmonly' type for background method ",background.method)
    default.b.param["type"] <- "together"
    background.param["type"]  <- "together"
  }

  for (name in names(background.param)){
    default.b.param[name] <- background.param[name]
  }

  
  if (is.element(as.character(response.variable),c("mm","MM"))& is.element(background.method,c("IdealMM","MASIM"))){
    if (R.model$mmorpm.covariate !=0){
      stop("Can't use a covariate PM with an Ideal Mismatch background")
    }
    warning("Ideal Mismatch correction will treat PM as MM and MM as PM")
    default.b.param["ideal"] <- "PM"
  } else if (is.element(as.character(response.variable),c("pm","PM"))& is.element(background.method,c("IdealMM","MASIM"))){
    if (R.model$mmorpm.covariate !=0){
      stop("Can't use a covariate MM with an Ideal Mismatch background")
    }
    default.b.param["ideal"] <- "MM"
  }
  

  
  default.b.param

}


verify.norm.param <- function(R.model, normalize.method,normalize.param = list()){


  get.default.parameters <- function(normalize.method){

    if (normalize.method == "quantile"){
      default.n.param <- list(type="separate")
    } else if (normalize.method == "quantile.probeset"){
      default.n.param <- list(type="separate",use.median=FALSE,use.log2=TRUE)
    } else if (normalize.method == "scaling"){
      default.n.param <- list(type="separate",scaling.baseline=-4,scaling.trim=0.0,log.scalefactors=FALSE)
    } else if (normalize.method == "quantile.robust"){
      default.n.param <- list(type="separate",use.median=FALSE,use.log2=FALSE,weights=NULL,remove.extreme = "variance", n.remove = as.integer(1))
    }    
  }


  validate.supplied.parameters <- function(normalize.method,supplied.parameters,defaults){
    if (!all(is.element(names(supplied.parameters),names(defaults)))){
      stop("At least one of the supplied normalization parameters is not known for this normalization method")
    }

    defaults[names(supplied.parameters)] <- supplied.parameters
    if (!is.element(defaults["type"],c("separate","pmonly","mmonly","together"))){
      stop("Supplied option",defaults["type"]," for 'type' is not valid")
    }


      
    if (normalize.method == "quantile.probeset"){
      if (!is.logical(defaults[["use.median"]])){
        stop("use.median should be TRUE or FALSE")
      }
      if (!is.logical(defaults[["use.log2"]])){
        stop("use.log2 should be TRUE or FALSE")
      }
    } else if (normalize.method == "scaling"){
      if (defaults[["scaling.trim"]] < 0 || defaults[["scaling.trim"]] >= 0.5){
        stop("scaling.trim can't be less than 0 or above 0.5")
      }
      if (defaults[["scaling.baseline"]] < -4){
        stop("scaling.baseline may be invalid")
      }

    } else if (normalize.method == "quantile.robust"){
      if (!is.logical(defaults[["use.median"]])){
        stop("use.median should be TRUE or FALSE")
      }
      if (!is.logical(defaults[["use.log2"]])){
        stop("use.log2 should be TRUE or FALSE")
      }
      if (!is.element(defaults[["remove.extreme"]],c("both","variance","mean"))){
        stop("remove.extreme ",defaults[["remove.extreme"]]," is not valid setting")
      }
      defaults["n.remove"] <- as.integer(defaults["n.remove"])

      if (is.character(defaults[["weights"]])){
        if (defaults[["weights"]] != "huber"){
          stop("The supplied weights option ",defaults[["weights"]]," is not valid")
        }
      } else if (is.double(defaults[["weights"]])){
        if (any(defaults[["weights"]] < 0)){
          stop("Can't have negative normalization weights")

        }
        if (sum(defaults[["weights"]] > 0) < 1) {
          stop("Need at least one non negative weights\n")
        }
      } else if (!is.null(defaults[["weights"]])){
        stop("Problem with the supplied weights option. Doesn't look valid.")
      }

      
    }
    defaults
  }


  
  

  oligoPLM.norm.methods <- c("quantile","scaling","quantile.probeset","quantile.robust")
  
  if (!is.element(normalize.method,oligoPLM.norm.methods)){
    stop(paste("Don't know the normalization method",normalize.method,"Please use one of the known methods:","quantile","scaling","quantile.probeset","quantile.robust",sep=" "))
  }

  default.n.param <- get.default.parameters(normalize.method)
  
  if (R.model$response.variable !=0){
    if (R.model$mmorpm.covariate == 0){
      if ( R.model$response.variable  == -1){
        default.n.param["type"] <- "mmonly"
      } else {
        default.n.param["type"] <- "pmonly"
      }
    }
  }

  if (R.model$response.variable == 0){
    if (is.element(default.n.param["type"],c("pmonly","mmonly"))){
      stop("Can't normalize 'pmonly' or 'mmonly' with PMMM response")
    }
  }

  default.n.param <- validate.supplied.parameters(normalize.method, normalize.param,  default.n.param) ####[names(normalize.param)] <- normalize.param

  if (is.element(normalize.param["type"],c("separate","together")) & R.model$response.variable !=0){
    if (R.model$mmorpm.covariate == 0){
      if ( R.model$response.variable  == -1){
        warning("Changing type in normalization to 'mmonly'")
        default.n.param["type"] <- "mmonly"
      } else {
        warning("Changing type in normalization to 'pmonly'")
        default.n.param["type"] <- "pmonly"
      }
    }
  }
    
  default.n.param
}



fitPLM <- function(object, model=PM ~ -1 + probes + samples,
                   variable.type=c(default="factor"),
                   constraint.type=c(default="contr.treatment"),
                   subset=NULL, background=TRUE, normalize=TRUE,
                   background.method="RMA.2",normalize.method="quantile",
                   background.param=list(), normalize.param=list(),
                   output.param=NULL,
                   model.param=NULL,
                   verbosity.level=0){

  if (!is(object, "FeatureSet")) {
    stop(paste("argument is", class(object), "fitPLM requires FeatureSet"))
  }

  pkgname <- annotation(object)
  if (!requireAnnotation(pkgname, verbose=TRUE))
      stop("The annotation package, ", pkgname, ", could not be loaded")

  if (is.null(output.param)) output.param <- verify.output.param()
  if (is.null(model.param)) model.param <- verify.model.param(object, model)
  
  b.param <- background.param
  n.param <- normalize.param
  
  variable.type <- verify.variable.types(model, variable.type)
  constraint.type <- verify.constraint.types(model, constraint.type)

  output <- verify.output.param(output.param)
  modelparam <- verify.model.param(object, model, model.param=model.param)
  R.model <- PLM.designmatrix3(object, model, variable.type=variable.type, constraint.type=constraint.type)

  background.param <- verify.bg.param(R.model, background.method, background.param = background.param)
  normalize.param <- verify.norm.param(R.model, normalize.method, normalize.param=normalize.param)
   
  if (!is.null(subset)){
    n.probesets <- length(subset)
  } else {
    n.probesets <- length(probesetNames(object))
  }

  ## to avoid having to pass location information to the c code, we will just call the R code method
  if (is.element(background.method, c("MAS", "MASIM")) & background){
    object <- backgroundCorrect(object, method="mas",
                                verbose=verbosity.level > 0)
  }
  if (is.element(background.method, c("gcrma", "GCRMA")) & background){
      stop("background correction via gcrma not yet implemented.")
  }
  pms <- pm(object, subset)

  if ('mmfeature' %in% dbListTables(db(get(pkgname)))){
      mms <- mm(object, subset)
  }else{
      mms <- matrix(0, 0, 0)
  }
  pns <- probeNames(object, subset)
  i <- order(pns)
  pms <- pms[i,, drop=FALSE]
  if (nrow(mms) != nrow(pms)){
      mms <- pms
  }else{
      mms <- mms[i,, drop=FALSE]
  }
  pns <- pns[i]
  rm(i)

  
  Fitresults <- .Call("R_rlm_PLMset_c",
                      pms, mms, pns,
                      length(unique(pns)),
                      R.model,
                      output,
                      modelparam,
                      background,
                      background.method,
                      background.param,
                      normalize,
                      normalize.method,
                      normalize.param,
                      verbosity.level,
                      PACKAGE="oligo")
  
  new("PLMset",
      chip.coefs=Fitresults[[1]],
      probe.coefs= Fitresults[[2]],
      weights=Fitresults[[3]],
      se.chip.coefs=Fitresults[[4]],
      se.probe.coefs=Fitresults[[5]],
      const.coefs=Fitresults[[6]],
      se.const.coefs=Fitresults[[7]],
      residuals=Fitresults[[8]],
      residualSE=Fitresults[[9]],
      varcov = Fitresults[[10]],
      cdfName = pkgname,
      phenoData = phenoData(object),
      annotation = pkgname,
      experimentData = experimentData(object),
      ##FIXME: remove # after notes is fixed.
      ##notes = object@notes,
      nrow= geometry(object)[1],
      ncol= geometry(object)[2],
      narrays=ncol(object),
      model.description = list(which.function="fitPLM",
        preprocessing=list(bg.method=background.method,bg.param=background.param,
          background=background,norm.method=normalize.method,
          norm.param=normalize.param,normalize=normalize),
        modelsettings =list(constraint.type=constraint.type,
          variable.type=variable.type,model.param=modelparam),
        outputsettings=output, R.model=R.model),
      manufacturer=tolower(manufacturer(object)))

}
