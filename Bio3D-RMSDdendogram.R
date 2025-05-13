library(bio3d)
#demo("pdb")
#demo("pca")

# Download some example PDB files
ids <- c("1XB7_A","2PJL_A","3D24_A","2PJL_B","3D24_C","3K6P_A","7E2E_A","7E2E_B")
raw.files <- get.pdb(ids)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

# Extract and align sequences
pdbs <- pdbaln(files)

core <- core.find(pdbs)

# See Figure 3.
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)

core.inds <- print(core, vol=1.0)

write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core.pdb")

xyz <- pdbfit( pdbs, core.inds )

rd <- rmsd(xyz)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd, k=3, labels=pdbs$id, cex=0.5,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
