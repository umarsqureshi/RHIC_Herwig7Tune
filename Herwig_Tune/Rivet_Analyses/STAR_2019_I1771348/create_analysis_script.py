
dirs = [f"{i:04d}" for i in range(100)]

out = open(f"generate_rivet_analyses.sh", 'w')

for dir in dirs:
    rivet_run = f"rivet --pwd -a STAR_2019_I1771348 --ignore-beams -o ../../out/{dir}/2019_{dir}.yoda  ../../out/{dir}/RHIC_Dijet.hepmc\n"
    out.writelines(rivet_run)
    
out.close()

