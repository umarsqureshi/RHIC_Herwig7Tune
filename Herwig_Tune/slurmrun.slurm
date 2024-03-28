#!/bin/bash

#SBATCH --account jswhep  # Which group to use for the job.
#SBATCH -J Herwig        # Job name
#SBATCH -o Herwig.o%j     # Name of stdout output file(%j expands to jobId)
#SBATCH -e Herwig.e%j     # Name of stderr output file(%j expands to jobId)
#SBATCH --mail-user=umar.sohail.qureshi@vanderbilt.edu
#SBATCH --constraint=zen
## mailtype begin, end, fail, all, time_limit_50, time_limit_80, time_limit_90
##SBATCH --mail-type=begin  # email me when the job starts
##SBATCH --mail-type=end    # email me when the job finishes
#SBATCH --mail-type=all    # email me when the job starts, ends, or fails
#SBATCH --mail-type=time_limit_80  # emal me when80% of time used
##SBATCH -p debug         # submit to the normal queue, partition production, debug, maxwell, pascal
#SBATCH --mem 32G         # amount of memory per node
#SBATCH --nodes=1        # Like -24, number of nodes on which to run
#SBATCH --ntasks=6
##SBATCH --ntasks-per-node=6
#SBATCH -t 12:40:00         # Run time (d-hh:mm:ss) 

echo "Preparing:"
date
set -x                          # Output commands
set -e                          # Abort on errors

source setup_herwig.sh

cd /home/sohailu/Herwig_Tune/out/0001
Herwig read RHIC_Diject.in
Herwig run RHIC_Diject.run -N 1000000






echo "Stopping:"
date
echo "Done."
