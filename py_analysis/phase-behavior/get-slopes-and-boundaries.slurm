#!/bin/bash
#
#SBATCH --job-name=boundary-maker-slope-finder
#SBATCH --ntasks=9                                     # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                              # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                              # total memory requested
#SBATCH --time=00:15:00                                # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                                # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu              # email address

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

for EMSA in -1 -0.5 0; do
	for EMSN in -1 -0.5 0; do
		srun --nodes=1 --ntasks=1 --mem-per-cpu=1GB --exclusive python boundary-probes.py --EMSA $EMSA --EMSN $EMSN > std_output_error_${EMSA}_${EMSN}.out 2>&1 &
	done
done

wait

