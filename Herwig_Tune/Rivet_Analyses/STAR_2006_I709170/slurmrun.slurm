#!/bin/bash

#SBATCH --account rke_group  # Which group to use for the job.
#SBATCH -J Rivet        # Job name
#SBATCH -o Rivet.o%j     # Name of stdout output file(%j expands to jobId)
#SBATCH -e Rivet.e%j     # Name of stderr output file(%j expands to jobId)
#SBATCH --mail-user=umar.sohail.qureshi@vanderbilt.edu
#SBATCH --constraint=zen
## mailtype begin, end, fail, all, time_limit_50, time_limit_80, time_limit_90
##SBATCH --mail-type=begin  # email me when the job starts
##SBATCH --mail-type=end    # email me when the job finishes
#SBATCH --mail-type=all    # email me when the job starts, ends, or fails
#SBATCH --mail-type=time_limit_80  # emal me when80% of time used
##SBATCH -p debug         # submit to the normal queue, partition production, debug, maxwell, pascal
#SBATCH --mem 8G         # amount of memory per node
#SBATCH --nodes=1        # Like -24, number of nodes on which to run
#SBATCH --ntasks=6
##SBATCH --ntasks-per-node=6
#SBATCH -t 1:40:00         # Run time (d-hh:mm:ss) 

echo "Preparing:"
date
set -x                          # Output commands
set -e                          # Abort on errors

cd

source setup_professor_with_rivet318.sh

cd /data/rke_group/Herwig_Tune/Rivet_Analyses/STAR_2006_I709170

./generate_rivet_analyses.sh


echo "Stopping:"
date
echo "Done."
