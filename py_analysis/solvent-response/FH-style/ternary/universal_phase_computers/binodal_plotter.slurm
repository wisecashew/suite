#!/bin/bash
#
#SBATCH --job-name=binodal-solver
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=00:30:00
#SBATCH --mail-type=all
#SBATCH --mail-user=satyend@princeton.edu

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
echo "Inside block 1..."

# srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive python bin_fc_v0.py --chisc 2.7 --chips 2.7 --chipc 2.7 -vc 1 -vs 1 -vp 1
# for c in 2.6 2.62 2.63 2.657 2.658 2.659 2.65 2.661 2.662 2.663 2.6666 2.667 2.668 2.669 2.66 2.68; do
for c in 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1; do
	srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive python compute.py --chisc $c --chips $c --chipc $c -vs 1 -vc 1 -vp 1 --plot-crits > output_${c}_1.out 2>&1 &
done

wait

echo "Number of computations is $x."
