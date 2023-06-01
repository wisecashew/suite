#!/bin/bash
#
#SBATCH --job-name=mu-solver
#SBATCH --ntasks=1        	     			                         # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                    			                 # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                            			         # total memory requested
#SBATCH --time=20:00:00                       			                 # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			                 # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			                 # email address

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

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
# python better_chem_pots.py --chiac 0 --chiab -1 --chibc -10 --mesh 300 -N 32
python find_binodals.py --dumpfile mu_bigger_sweep.txt --image binodals-extra > output
