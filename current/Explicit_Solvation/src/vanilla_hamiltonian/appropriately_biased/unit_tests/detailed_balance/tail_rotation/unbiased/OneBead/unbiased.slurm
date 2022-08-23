#!/bin/bash
#
#SBATCH --job-name=unit_test_unbiased_singlebead                           # name of job
#SBATCH --ntasks=12        	     			                   # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                    			           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                            			   # total memory requested
#SBATCH --time=00:30:00                       			           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      			           # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    			           # email address

# The modules I need to run my job
module purge
# module load rh/devtoolset/7 

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?." ' EXIT

for dir in T_2 T_3; do
    cd $dir
    for i in {1..3}; do
        ./hr_o_n_boltz -f 100000 -M 10'000'000 -p initial.txt -t geom_and_esurf.txt -u energy1_$i -L lattice1_$i -e orientation1_$i -s stats1_$i -o coords1_$i > o_to_n_$i 2>&1 &
        ./hr_n_o_boltz -f 100000 -M 10'000'000 -p final.txt -t geom_and_esurf.txt -u energy2_$i -L lattice2_$i -e orientation2_$i -s stats2_$i -o coords2_$i > n_to_o_$i 2>&1 & 
    done 
    cd ..
done 

echo "All simulations have been sent..."

wait 

echo "All simlations done!"
