#!/bin/bash
#
#SBATCH --job-name=N:32-entropy_regime-appropriately-biased                           # name of job
#SBATCH --ntasks=1        	     			                   # number of tasks across all nodes 
#SBATCH --cpus-per-task=48                    			           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=100MB                            			   # total memory requested
#SBATCH --time=00:05:00                       			           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			           # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			           # email address

export OMP_NUM_THREADS=48
# The modules I need to run my job
module purge

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~

rm output.mc energydump.mc
./cosolv_ani -f 1 -M 20 -p 32mer.txt -t geom_and_esurf.txt -o coords.mc -L lattice_dump.mc -u energydump.mc -e orientation.mc -s stats.mc > output.mc 2>&1

echo "running simulations..."

