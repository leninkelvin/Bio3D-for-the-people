library(bio3d)
library(httr)

# Download PDB, calcualte normal modes of the open subunit
pdb.full   <- read.pdb("2LAO")
pdb.other   <- read.pdb("1LST")
pdb.open   <- trim.pdb(pdb.full, atom.select(pdb.full, chain="A"))
modes      <- nma(pdb.open)

print(modes)

plot(modes, sse=pdb.open)

# Calculate modes with various force fields
modes.a <- nma(pdb.open, ff="calpha")
modes.b <- nma(pdb.open, ff="anm")
modes.c <- nma(pdb.open, ff="pfanm")
modes.d <- nma(pdb.open, ff="reach")
modes.e <- nma(pdb.open, ff="sdenm")

# Root mean square inner product (RMSIP)
r <- rmsip(modes.a, modes.b)

plot(r, xlab="ANM", ylab="C-alpha FF")

# Make a PDB trajectory
mktrj(modes, mode=7)

# Vector field representation (see Figure 3.)
pymol(modes, mode=7)

# Calculate the cross-correlation matrix
cm <- dccm(modes)

# Plot a correlation map with plot.dccm(cm)
plot(cm, sse=pdb.open, contour=F, col.regions=bwr.colors(20), at=seq(-1,1,0.1) )

# View the correlations in the structure (see Figure 5.)
#pymol(cm, pdb.open, type="launch")

# Deformation energies
defe <- deformation.nma(modes)
defsums <- rowSums(defe$ei[,1:3])

# Fluctuations
flucts <- fluct.nma(modes, mode.inds=seq(7,9))

# Write to PDB files (see Figure 6.)
write.pdb(pdb=NULL, xyz=modes$xyz, file="R-defor.pdb", b=defsums)
write.pdb(pdb=NULL, xyz=modes$xyz, file="R-fluct.pdb", b=flucts)

# Closed state of the subunit
pdb.closed <- trim.pdb(pdb.full, atom.select(pdb.other, chain="A"))

# Align closed and open PDBs
aln <- struct.aln(pdb.open, pdb.closed, max.cycles=0)
pdb.closed$xyz <- aln$xyz

# Caclulate a difference vector
xyz <- rbind(pdb.open$xyz[aln$a.inds$xyz], pdb.closed$xyz[aln$a.inds$xyz])
diff <- difference.vector(xyz)

# Calculate overlap
oa <- overlap(modes, diff)

plot(oa$overlap, type='h', xlab="Mode index", ylab="Squared overlap", ylim=c(0,1))
points(oa$overlap, col=1)
lines(oa$overlap.cum, type='b', col=2, cex=0.5)
text(c(1,5)+.5, oa$overlap[c(1,5)], c("Mode 1", "Mode 5"), adj=0)

# Run geostas to find domains
gs <- geostas(modes, k=4)

# Write NMA trajectory with domain assignment
mktrj(modes, mode=7, chain=gs$grps)

# Define the ensemble PDB-ids
ids <- c("1XB7_[A]", "2PJL_[A]", "3D24_[A]", "3K6P_[A]", "7E2E_[A]")

# Download and split PDBs by chain ID
raw.files <- get.pdb(ids, "groel_pdbs", gzip=TRUE)
files <- pdbsplit(raw.files, ids, path = "groel_pdbs")

# Align and superimpose coordinates
pdbs <- pdbaln(files, fit=TRUE)

# Run geostast to find domains
gs <- geostas(pdbs, k=4)

# Plot a atomic movement similarity matrix
plot(gs, margin.segments=gs$grps, contour=FALSE)

# Principal component analysis 
gaps.pos <- gap.inspect(pdbs$xyz)
xyz <- fit.xyz(pdbs$xyz[1, gaps.pos$f.inds],
               pdbs$xyz[, gaps.pos$f.inds],
               fixed.inds=gs$fit.inds,
               mobile.inds=gs$fit.inds)

pc.xray <- pca.xyz(xyz)

# Visualize PCs with colored domains (chain ID)
mktrj(pc.xray, pc=1, chain=gs$grps)
