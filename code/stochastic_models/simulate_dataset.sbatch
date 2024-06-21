#!/bin/bash

module load R/3.5.3
#Rscript "/michorlab/jdean/p21_central_dogma/code/analyses/simulate_dataset_p21rna.R" "${arg1}" "${arg2}" "${arg3}" "${arg4}" "${arg5}"
Rscript "/michorlab/jdean/p21_central_dogma/code/stochastic_models/simulate_dataset_p21rna.R" "${arg1}" "${arg2}" "${arg3}" "${arg4}" "${arg5}" "${arg6}" "${arg7}" "${arg8}" "${arg9}" "${arg10}" "${arg11}"
