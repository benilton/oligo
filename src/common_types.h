/*********************************************************************
 **
 ** file: common_types.h
 **
 ** Aim: Define some structures that will be used for passing
 **      related data items to model fitting routines.
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Sept 04, 2003
 **
 ** History
 ** Sept 04, 2003 - Initial version (outputsettings, DataGroup, 
 **                 PLMoutput, RPLMoutput, PLMmodelparam)      
 ** Sept 13, 2003 - Added n_rlm_iterations and init_method to
 **                 PLMmodelparam (these will control how many
 **                 iterations of iteratively reweighted least
 **                 squares and how we intialize the IRLS: least 
 **                 squares or median polish.
 ** May 27, 2004  - added default.model to  PLMmodelparam if non zero then
 **                 we are fitting -1 + samples + probes with (sum probes = 0 constraint)
 ** Feb 28, 2006 - change // to ansi c comments for old compilers
 **
 *********************************************************************/

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H 1

#include "threestep_summary_methods.h"




#include <R.h>
#include <Rinternals.h>

/*******************************************************************
 **
 ** a structure for holding flags and settings on what should
 ** be stored as output from the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  int weights;   /* Store Weights */
  int residuals;     /* Store Residuals */
  int residSE;   /* Store ResidSE */
  int pseudoSE;  /* Compute Pseudo SE in the case of an rmaPLM */
  int varcov;   /* Store varcov matrices what type. 0 means none, 1 means only chip level, 2 means all,*/
} outputsettings;



/*******************************************************************
 **
 ** a structure for holding probe intensity data as it gets
 ** used by the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  double *PM;
  double *MM;
  int rows;
  int cols;
  int nprobesets;
  const char **ProbeNames;
} Datagroup;

/*******************************************************************
 **
 ** a structure for holding computed quantities that result from
 ** using the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  char **outnames;
  double *out_weights;
  double *out_probeparams;
  double *out_chipparams;
  double *out_constparams;
  double *out_probe_SE;
  double *out_chip_SE;
  double *out_const_SE;
  double *out_resids;
  double *out_residSE;
  double **out_varcov;

} PLMoutput;


/*******************************************************************
 **
 ** a structure for holding computed quantities as they will be 
 ** returned to R (ie in R data structures
 **
 **
 **
 ******************************************************************/

typedef struct {
  SEXP weights;
  SEXP probe_coef; 
  SEXP chip_coef;
  SEXP const_coef;
  SEXP chip_SE;
  SEXP probe_SE;
  SEXP const_SE;
  SEXP residuals;
  SEXP residSE;
  SEXP varcov;
  int nprotected;
} PLMRoutput;

/*******************************************************************
 **
 ** a structure for holding parameters and flags that describe the
 ** PLM model being fitted.
 ** 
 **
 **
 **
 ******************************************************************/

typedef struct{
  int nchipparams;
  int method;  
  int se_method;
  int psi_code;
  double psi_k;
  double *input_chipcovariates;
  int n_rlm_iterations;
  int init_method;
  int default_model;
  int mmorpm_covariate; 
  pt2PLMSummary PLM3stepSummary;
} PLMmodelparam;




/******************************************************************
 ** Below here are new structures
 ******************************************************************/



/*******************************************************************
 **
 ** a structure for holding parameters and flags that describe the
 ** PLM model being fitted. These are universal irrespective of 
 ** the current probeset (ie the number of probes does not change
 ** 
 **
 **
 ******************************************************************/

typedef struct{
  int psi_code;   /* A specifier for which M-estimator to use */
  int se_method;  /* A specified for which SE method to use */
  double psi_k;   /* A parameter used for control of M-estimator influence function */
  int n_rlm_iterations; /* maximum number of iteratively reweighted least squares iterations */
  int init_method; /* should we start IRLS with least squares or Huber  */
  int mmorpm_covariate; /* -1 means MM is response PM is covariate, 0 Means no PM or MM covariate, 1 means PM is response MM is covariate */
  int response_variable; /* -1 means MM response, 0 means both PM and MM are both response variables,  1 means PM is the response */
  int *which_parameter_types;   /* a vector of length 5 a zero element in an element means that that that parameter is not in model 
				** element 0   - intercept?
				** element 1   - is there is chip-level treatment factor and covariate variables?
				** element 2   - samples (ie chip effect)
				** element 3   - probe.types effect (note this is only valid for PMMM response variables)
				** element 4   - probe effect */
  int *strata;                  /* a vector of length 5, each element corresponds to the parameter order above
			        ** note that only probe.types and probe.effect are computed within strata */
  int *constraints;             /* a vector of length 5 specifying the constraint types for parameters same order as specified above
				**  -1 = sum to zero 
                                **   0 = unconstrained
                                **   1 = first treatment = 0 */
  int *probe_type_treatment_factor;  /* a vector of the same length as the number of arrays. values should be between 0, and max_probe_type_treatment_factor
				     ** used if the probe type effect is to be estimated within a treatment factor strata  */
  int max_probe_type_treatment_factor;
  int *probe_treatment_factor;       /* a vector of the same length as the number of arrays. values should be between 0 and max_probe_treatment_factor
				     ** used if the probe effect is to be estimated within  a treatment factor strata */
  int max_probe_treatment_factor;
  double *chiplevelcovariates;       /* a vector of values that are covariates (or treatment/genotype factor) variables related to chips */
  int n_chiplevelcovariates;
  int n_arrays;                      /* the number of arrays */
  double *input_chip_weights;        /* NULL if no input weights, otherwise length is n_arrays and values are non negative, should be scaled to be on 0-1 */
  double *input_probe_weights;       /* NULL if no input weights, otherwise length is data->n_probes or 2*data->n_probes (depending on if PM,MM or PMMM model)
				     ** should be scaled to be on 0-1 and non-negative */
  int trans_fn;                      /* a code which corresponds to how data should be transformed */
} PLM_model_parameters;


/*******************************************************************
 **
 ** a structure for holding probe intensity data as it gets
 ** used by the PLM routines.
 **
 **
 **
 ******************************************************************/

typedef struct {
  double *PM;           /* a matrix of perfect match intensities   dim:nprobes by n_arrays */
  double *MM;           /* a matrix of mismatch intensities        dim:nprobes by n_arrays */
  int n_probesets;       /* the number of probesets */
  int n_arrays;          /* the number of arrays */
  int n_probes;          /* number of probes on each array */
  const char **ProbeNames; /* name of the probeset to which the probeset belongs */
} PLM_Datagroup;


/**********************************************************************
 **
 ** the modelfit struct is used for storing information about the 
 ** current model (the one being fitted individually to each probeset.
 ** It stores the model matrix and places to output the model
 **
 **********************************************************************/

typedef struct{
  double *cur_params;            /* storage for output */
  double *cur_se_estimates;
  double *cur_weights;
  double *cur_resids;
  double *cur_varcov;
  double *cur_residSE;
  double *X;      /* design matrix */
  int n;          /* number of observations */
  int p;          /* number of parameters */
  int nprobes;    /* number of probes in current probeset */
} PLM_modelfit;


/**********************************************************************
 **
 ** the output struct is used for storing all the output from
 ** the model fitting procedure
 **
 **
 **********************************************************************/

typedef struct {
  char **outnames;   /* Names of each probeset */



  /* the following is where estimated parameters are stored */
  double **out_probeparams; /* all probe effect parameters */
  double *out_chipparams;  /* Chip level factor/covariate and sample-effect values */
  double *out_constparams; /* intercept, MM/PM covariate, probe.types */

  /* the following is where the estimated SE are stored */
  double **out_probe_SE;
  double *out_chip_SE;
  double *out_const_SE; 

  
  double **out_weights; /* weights is two matrices the first matrix is PM
			** the second matrix is MM */
  double **out_resids;  /* residuals are two matrices the first matrix is PM
			** the second matrix is MM */
  double *out_residSE;
  double **out_varcov;

} PLM_output;


/**********************************************************************
 **
 ** this struct specifies whether the optional components should
 ** be stored
 **
 **
 **
 *********************************************************************/

typedef struct {
  int weights;   /* Store Weights */
  int residuals;     /* Store Residuals */
  int residSE;   /* Store ResidSE */
  int varcov;   /* Store varcov matrices what type. 0 means none, 1 means only chip level, 2 means all,*/
} PLM_outputsettings;






PLM_modelfit *new_PLMmodelfit();
void free_PLMmodelfit(PLM_modelfit *x);



#endif
