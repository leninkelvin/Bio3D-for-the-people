library(bio3d)
library(httr)

# Select bacterial DHFR PDB IDs
ids <- c("1XB7_A","2PJL_A","3D24_A","3K6P_A","7E2E_A")
#ids <- c("3WX1_A", "3WX2_A", "4TZC_A", 
#         "4TZU_A", "5YIZ_A", "5YJ0_A",
#         "5YJ1_A")

# Download and split by chain ID
raw.files <- get.pdb(ids, path="raw_pdbs")
files     <- pdbsplit( raw.files, ids )

# Alignment of structures
pdbs <- pdbaln(files)

# Sequence identity
summary( c(seqidentity(pdbs)) )

# NMA on all structures
modes <- nma(pdbs)

print(modes)

# Plot fluctuation data
col <- c(1, 2, 3, 4, 5)
plot(modes, pdbs=pdbs, col=col)

#legend("topleft", col=unique(col), lty=1,
       #legend=c("E.Coli", "B.Anthracis", "M.Tubercolosis", "S.Aureus"))

# Alternatively, one can use 'rm.gaps=FALSE' to keep the gap containing columns
modes <- nma.pdbs(pdbs, rm.gaps=FALSE)

# Calculate correlation matrices for each structure
cij <- dccm(modes)

# Determine correlations present only in all 5 input structures
cij.all <- filter.dccm(cij$all.dccm, cutoff.sims=5, cutoff.cij = 0)
plot.dccm(cij.all, main="Consensus Residue Cross Correlation")
