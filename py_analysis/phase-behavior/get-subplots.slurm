#!/bin/bash
#
#SBATCH --job-name=subplot-arguments
#SBATCH --ntasks=2        	     			                   # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                    			           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=3GB                            			   # total memory requested
#SBATCH --time=00:15:00                       			           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			           # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			           # email address

# The modules I need to run my job
module purge
module load anaconda3/2020.7
conda activate data_analysis

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~

srun --nodes=1 --ntasks=1 --mem-per-cpu=3GB --exclusive python heatmaps-vectorized-v2.py > heat-v2.out 2>&1 &
srun --nodes=1 --ntasks=1 --mem-per-cpu=3GB --exclusive python heatmaps-vectorized-longer-v2.py > heat-v2-longer.out 2>&1 &

wait

