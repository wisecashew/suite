#!/bin/bash
#
#SBATCH --job-name=binodal-solver
#SBATCH --ntasks=13
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
for c in 2.4 2.5 2.6 2.7 2.8 2.9 3.0; do
	srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive python compute.py --chisc $c --chips $c --chipc $c -vs 1 -vc 1 -vp 1 --plot-crits > output_$c.out 2>&1 &
done

for c in 8 9 10; do
	srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive python compute.py --chisc -$c --chips -1 --chipc -1 -vs 1 -vc 1 -vp 32 --plot-crits > output_$c.out 2>&1 &
done

for c in 8 9 10; do
	srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive python compute.py --chisc 0 --chips -1 --chipc -$c -vs 1 -vc 1 -vp 32 --plot-crits > output__$c.out 2>&1 &

done
wait

echo "Number of computations is $x."
