import os

dirList = os.listdir("./")
dirs = [dir for dir in dirList if dir.startswith("00")]


for dir in dirs:

    slurm_script = f"""#!/bin/bash

#SBATCH --account rke_group  # Which group to use for the job.
#SBATCH -J Her{dir}        # Job name
#SBATCH -o Herwig{dir}.o%j     # Name of stdout output file(%j expands to jobId)
#SBATCH -e Herwig{dir}.e%j     # Name of stderr output file(%j expands to jobId)
#SBATCH --mail-user=umar.sohail.qureshi@vanderbilt.edu
##SBATCH --constraint=zen
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
#SBATCH -t 12:00:00         # Run time (d-hh:mm:ss) 

echo "Preparing:"
date
set -x                          # Output commands
set -e                          # Abort on errors

cd

source setup_herwig.sh

cd /data/rke_group/Herwig_Tune/out/{dir}
Herwig read RHIC_Dijet.in
Herwig run RHIC_Dijet.run -N 10000000

echo "Stopping:"
date
echo "Done."
"""

    out = open(f"{dir}/slurmrun.slurm", 'w')
    out.writelines(slurm_script)
    out.close()