#!/bin/bash
#
#SBATCH --job-name=binodal-tracer
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=02:30:00
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
srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive --mem-per-cpu=2G python universal.py --chisc -9 --chips -1 --chipc -1 -vs 1 -vc 1 -vp 5 --final-binodal BINODALS/chips_-1-chipc_-1-chisc_-9-vs_1-vc_1-vp_5.binodals.pkl --search-density 1000 --island-stable-pkl STABLE_ISLANDS/chips_-1-chipc_-1-chisc_-9-vs_1-vc_1-vp_5.stable_islands.pkl --island-unstable-pkl UNSTABLE_ISLANDS/chips_-1-chipc_-1-chisc_-9-vs_1-vc_1-vp_5.unstable_islands.pkl --mesh-pkl MESHES/chips_-9-chipc_-1-chisc_-1-vs_1-vc_1-vp_1.mesh.pkl --crit-pkl CRITS/chips_-1-chipc_-1-chisc_-9-vs_1-vc_1-vp_5.crits.pkl --plot-binodals --plot-crits

echo "Number of computations is $x."
