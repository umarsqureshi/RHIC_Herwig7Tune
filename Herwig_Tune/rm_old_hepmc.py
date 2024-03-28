dirs = [f"{i:04d}" for i in range(100)]

out = open(f"remove_old_hepmc.sh", 'w')

for dir in dirs:
    rm = f"cd /home/sohailu/RHIC_Herwig7Tune/Herwig_Tune/out/{dir}\nrm -r RHIC_Dijet.hepmc\n"
    out.writelines(rm)
    
out.close()

