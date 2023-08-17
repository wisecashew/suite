#!/bin/bash


set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
declare -A mesh
mesh[7]=200
mesh[8]=225
mesh[9]=250
mesh[10]=275
mesh[11]=300
mesh[12]=325

# for chiac in {7..11}; do
# 	echo "In chiac = -${chiac}..."
# 	python pskelbin.py --chiac -${chiac} --chiab -1 --chibc -1 --mesh ${mesh[${chiac}]} -N 32 --skelfile binskel-N_32-chiab_-1-chibc_-1-chiac_-${chiac}.skelfile > chiac_-${chiac}.out 2>&1
# done

mesh[8]=200
mesh[9]=225
mesh[10]=250
mesh[11]=250
mesh[12]=250

for chibc in {11..12}; do
	echo "In chibc = -${chibc}..."
	python pskelbin.py --chiac 0 --chiab -1 --chibc -${chibc} --mesh ${mesh[${chibc}]} -N 32 --skelfile binskel-N_32-chiab_-1-chibc_-${chibc}-chiac_0.skelfile > chibc_-${chibc}.out 2>&1 
done

sbatch binodal_maker.slurm
