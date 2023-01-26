#!/bin/bash
#
#SBATCH --job-name=histogram-test-biased                           # name of job
#SBATCH --ntasks=6        	     			                   # number of tasks across all nodes 
#SBATCH --cpus-per-task=1                    			           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1GB                            			   # total memory requested
#SBATCH --time=06:30:00                       			           # total run time limit (HH:MM:SS)
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

for dir in biased; do
    cd $dir
    for i in {1..3}; do
        cd T_$i 
        echo "Inside $dir/T_$i"
        biased -f 1000 -M 10'000'000 -p 8mer.txt -t geom_and_esurf.txt -u energydump -L lattice -e orientation -s stats -o coords > output 2>&1 & 
        cd ..
    done 
    cd ..
done 

echo "All biased simulations have been sent..."

wait 

echo "All simlations done!"
