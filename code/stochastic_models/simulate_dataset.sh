#!/bin/bash

# location for home directory to work with
folder_dir="michorlab/jdean/p21_central_dogma"
# Change read/write privileges to allow me to edit my parameters file
# Each parameter set I run is a line
chmod 777 "$stochastic_models/parameters.dat"

while IFS=$',' read -r -a params
do
 # Command to send to nodes - will use qsub here
 # -v argument_name=argument_value will pass these arguments into run_multi_jobs3.sh

  #echo sbatch -o output.txt -e error.txt --job-name=myjob --mem=10000 --export=arg1=${params[0]},arg2=${params[1]},arg3=${params[2]},arg4=${params[3]},arg5=${params[4]} p21rna.sbatch
  #sbatch -o output.txt -e error.txt --job-name=myjob --mem=10000 --export=arg1=${params[0]},arg2=${params[1]},arg3=${params[2]},arg4=${params[3]},arg5=${params[4]} p21rna.sbatch

  sbatch -o output.txt -e error.txt --job-name=myjob --mem=10000 --export=arg1=${params[0]},arg2=${params[1]},arg3=${params[2]},arg4=${params[3]},arg5=${params[4]},arg6=${params[5]},arg7=${params[6]},arg8=${params[7]},arg9=${params[8]},arg10=${params[9]},arg11=${params[10]} simulate_dataset.sbatch

done < parameters.dat

