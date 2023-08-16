#!/bin/bash
#
#SBATCH --job-name=mu-solver
#SBATCH --ntasks=1        	     			                         # number of tasks across all nodes 
#SBATCH --cpus-per-task=96                    			                 # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                            			         # total memory requested
#SBATCH --time=02:30:00                       			                 # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			                 # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			                 # email address

# The modules I need to run my job
module purge
module load anaconda3/2020.7
conda activate phase
# conda activate data_analysis
# module load rh/devtoolset/7 

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~

srun --ntasks=1 --cpus-per-task=1 --exclusive --nodes=1 python spinwbin_v3.py --chiac -7 --chiab -1 --chibc -1 --nproc 96 -N 32 --dumpfile binskel-chiab_-1-chibc_-1-chiac_-7.skelfile --bin-boundary bin-chiab_-1-chibc_-1-chiac_-7.out --image binodal-chiab_-1-chibc_-1-chiac_-7.png
