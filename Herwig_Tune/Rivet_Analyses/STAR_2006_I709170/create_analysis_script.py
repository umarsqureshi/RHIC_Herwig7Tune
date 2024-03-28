
dirs = [f"{i:04d}" for i in range(100)]

out = open(f"generate_rivet_analyses.sh", 'w')

for dir in dirs:
    rivet_run = f"rivet --pwd -a STAR_2006_I709170 --ignore-beams -o ../../out/{dir}/pispectra{dir}.yoda  ../../out/{dir}/RHIC_Dijet.hepmc\n"
    out.writelines(rivet_run)
    
out.close()

