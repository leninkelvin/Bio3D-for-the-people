library(bio3d)
demo("pdb")
demo("pca")
library(muscle)
library(httr)

## WORKING WITH A SINGLE STRUCTURE

pdb <- read.pdb("1Z05")

pdb

ca.inds <- atom.select(pdb, "calpha")

# See Figure 1
bf <- pdb$atom$b[ca.inds$atom]
plot.bio3d(bf, resno=pdb, sse=pdb, ylab="B-factor", xlab="Residue", typ="l")

##WORKING WITH MULTIPLE PDBs

# Download some example PDB files
ids <- c("1XB7_A","2PJL_A","3D24_A","3K6P_A","7E2E_A")
raw.files <- get.pdb(ids)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

## Calculate RMSD
rmsd(pdbs, fit=TRUE)

#FINDING SIMILAR STRUCTURES

pdb <- read.pdb("1Z05")
seq <- pdbseq(pdb)
blast <- blast.pdb(seq)

# See Figure 2.
hits <- plot.blast(blast, cutoff=200)

head(hits$hits)

head(hits$pdb.id)

anno <- pdb.annotate(hits$pdb.id)
head(anno[, c("resolution", "ligandId", "citation")])

      annotation <- pdb.annotate(hits$pdb.id)
      head(annotation[, c("resolution", "ligandId", "citation")])

#MULTIPLE SEQUENCE ALIGNMENT
# Download PDBs and split by chain ID
files <- get.pdb(hits, path="raw_pdbs", split = TRUE)

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

#RMSF
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# Tailor the PDB structure to exclude gap positions for SSE annotation
id = grep("1Z05", pdbs$id)
ref.pdb = trim.pdb(pdb, inds=atom.select(pdb, resno = pdbs$resno[id, gaps.res$f.inds]))

# Plot RMSF with SSE annotation and labeled with residue numbers (Figure 8.)
rf <- rmsf(xyz[, gaps.pos$f.inds])
plot.bio3d(rf, resno=ref.pdb, sse=ref.pdb, ylab="RMSF (Å)", xlab="Residue No.", typ="l")

#tor <- torsion.pdb(pdb)

# Basic Ramachandran plot (Figure 9)
#plot(tor$phi, tor$psi, xlab="phi", ylab="psi")

#COMPARING TWO STRUCTURES
# Locate the two structures in pdbs
#ind.a <- grep("1XB7_A", pdbs$id)
#ind.b <- grep("7E2E_A", pdbs$id)

# Exclude gaps in the two structures to make them comparable
#gaps.xyz2 <- gap.inspect(pdbs$xyz[c(ind.a, ind.b), ])
#a.xyz <- pdbs$xyz[ind.a, gaps.xyz2$f.inds]
#b.xyz <- pdbs$xyz[ind.b, gaps.xyz2$f.inds]

# Compare CA based pseudo-torsion angles between the two structures
#a <- torsion.xyz(a.xyz, atm.inc=1)
#b <- torsion.xyz(b.xyz, atm.inc=1)
#d.ab <- wrap.tor(a-b)
#d.ab[is.na(d.ab)] <- 0

# Plot results with SSE annotation
#plot.bio3d(abs(d.ab), resno=pdb, sse=pdb, typ="h", xlab="Residue No.", ylab = "Difference Angle")

#DISTANCE MATRIX

#a <- dm.xyz(a.xyz)
#b <- dm.xyz(b.xyz)

#plot.dmat( (a - b), nlevels=10, grid.col="gray", xlab="1tag", ylab="1tnd")

# Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
pc.xray

plot(pc.xray)

plot(pc.xray, pc.axes=1:2)

# Left-click on a point to label and right-click to end
#identify(pc.xray$z[,1:2], labels=basename.pdb(pdbs$id))

par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC1")
plot.bio3d(pc.xray$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC2")
plot.bio3d(pc.xray$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC3")

# See Figure 15
mktrj.pca(pc.xray, pc=1, file="pc1.pdb")

hc <- hclust(dist(pc.xray$z[,1:2]))
grps <- cutree(hc, h=30)
cols <- c("red", "green", "blue")[grps]
plot(pc.xray, pc.axes=1:2, col=cols)

# Dendrogram plot
names(cols) <- pdbs$id
hclustplot(hc, colors=cols, ylab="Distance in PC Space", main="PC1-2", cex=0.5, fillbox=FALSE)

#identify(pc.xray$z[,1], pc.xray$z[,2], labels=pdbs$id)



#plot(ENERGY$V15, xlim = c(0,54)) # For MMGBSA
