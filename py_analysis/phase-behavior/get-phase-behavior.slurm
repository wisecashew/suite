#!/bin/bash
#
#SBATCH --job-name=phase-makers
#SBATCH --ntasks=2                                     # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                              # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=100GB                              # total memory requested
#SBATCH --time=02:00:00                                # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                                # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu              # email address

# The modules I need to run my job
module purge
module load anaconda3/2020.7
conda activate phase

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
### code block for running simulations in the simple FH regime 
#~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~

srun --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=100GB python h4hourglass.py --png-name spinodal_corr_tot -pv 1.0 --delta-ms 0 -0.523 -0.5237 -1.3 > outputcorr0.net 2>&1 &
srun --nodes=1 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=100GB python h4hourglass.py --png-name spinodal_lock_tot -pv 0.25 --delta-ms -1 -2.07 -2.075 -3.7 > outputlock0.net 2>&1 &

wait
# python h4hourglass.py --png-name spinplots_net  --range-index 0 > output0.net 2>&1 &
# python h4hourglass.py --png-name spinplots_neck --range-index 1 > output1.net 2>&1 &

# wait

# for i in 0.0 0.1 0.25 0.5 1.0; do 
# 	python phase-dw.py --pv $i --png-name phases-pv-$i > output_$i.phase 2>&1 &
# done

# wait
