#!/bin/bash
#
#SBATCH --job-name=MCPolymerIsing,kT=1/2            # create a short name for your job
#SBATCH --qos=short              	     # _quality of service_
#SBATCH --nodes=1                     	     # node count
#SBATCH --ntasks-per-node=1          	     # number of tasks per node
#SBATCH --cpus-per-task=1                    # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4GB                            # total memory requested
##SBATCH --mem-per-cpu=4G                    # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=0:50:00                       # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all                      # send email on job start, end, and fail
#SBATCH --mail-user=satyend@princeton.edu    # email address
#
# The modules I need to run my job
module purge
module load rh/devtoolset/7 

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# echo an error message before exiting
# trap 'echo "\"${last_command\" command filed with exit code $?." ' EXIT
rm coords.txt energydump.txt stderr.txt

echo "got rid of the possibly problematic files?" 

echo "executing monte carlo..." 

./MonteCarlo -f 100 -N 100000 -p positions_8.txt -t geom_and_esurf.txt -o coords.txt -v > stderr.txt 2>&1  

echo "performed the execution of MC Polymer Ising!" 

