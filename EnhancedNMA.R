library(bio3d)

# read PDB
pdb <- read.pdb("1XB7")

# keep only protein + methotrexate
pdb <- trim(pdb, "notwater")

# calculate all-atom NMA with ENM, output calpha only
m.aa <- aanma(pdb, outmodes="calpha")
# write summary
m.aa
# compare with c-alpha modes
m.ca <- nma(pdb)
rmsip(m.aa, m.ca)
# all-atom NMA, but output all atoms
m.noh <- aanma(pdb, outmodes="noh")
# output trajectory of first modes
mktrj(m.noh, mode = 7, pdb=pdb)
