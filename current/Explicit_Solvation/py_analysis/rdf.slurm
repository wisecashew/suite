#!/bin/bash
#
#SBATCH --job-name=rdfs                          # name of job
#SBATCH --ntasks=1        	     			                   # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                    			           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                            			   # total memory requested
#SBATCH --time=00:15:00                       			           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			           # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			           # email address

# The modules I need to run my job
module purge
module load anaconda3/2020.7
conda activate data_analysis
# module load rh/devtoolset/7 

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

python plot_rdf_s1s2.py > rdf-plotter 2>&1 

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~


echo "Executed all!" 
echo "Executed DOP = 32 for all!" >> README
echo "Slurm job id is ${SLURM_JOB_ID}." >> README

