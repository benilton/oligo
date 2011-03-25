#ifndef PLM_MODELMATRIX_H
#define  PLM_MODELMATRIX_H 1

#include "common_types.h"

void PLM_build_model_matrix(const PLM_model_parameters *model, const PLM_Datagroup *data, PLM_modelfit *current, int *current_rows, int new_nprobes);

#endif
