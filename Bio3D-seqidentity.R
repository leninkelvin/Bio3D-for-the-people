library(bio3d)
#demo("pdb")
#demo("pca")

# Download some example PDB files
#ids <- c("3WX1_A", "3WX2_A", "4TZC_A", 
#         "4TZU_A", "5YIZ_A", "5YJ0_A",
#
ids <- c("1XB7_A","2PJL_A","3D24_A","2PJL_B","3D24_C","3K6P_A","7E2E_A","7E2E_B")
raw.files <- get.pdb(ids)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

