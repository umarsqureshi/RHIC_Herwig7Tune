dirs = [f"{i:04d}" for i in range(100)]

out = open(f"remove_old.sh", 'w')

for dir in dirs:
    rm = f"cd /data/rke_group/Herwig_Tune/out/{dir}\nrm -r STAR_2006_I709170.yoda\n"
    out.writelines(rm)
    
out.close()

