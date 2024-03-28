import os

dirList = os.listdir("./")
dirs = [dir for dir in dirList if dir.startswith("00")]

out = open(f"generate_MC.sh", 'w')

for dir in dirs:

    slurm_script = f"""cd /data/rke_group/Herwig_Tune/out/{dir}
sbatch slurmrun.slurm
"""
    out.writelines(slurm_script)
    
out.close()