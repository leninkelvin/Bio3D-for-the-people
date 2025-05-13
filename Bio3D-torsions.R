library(bio3d)

#COMPARING TWO STRUCTURES
# Locate the two structures in pdbs
ind.a <- grep("1XB7_A", pdbs$id)
ind.b <- grep("2PJL_A", pdbs$id)

# Exclude gaps in the two structures to make them comparable
gaps.xyz2 <- gap.inspect(pdbs$xyz[c(ind.a, ind.b), ])
a.xyz <- pdbs$xyz[ind.a, gaps.xyz2$f.inds]
b.xyz <- pdbs$xyz[ind.b, gaps.xyz2$f.inds]

# Compare CA based pseudo-torsion angles between the two structures
a <- torsion.xyz(a.xyz, atm.inc=1)
b <- torsion.xyz(b.xyz, atm.inc=1)
d.ab <- wrap.tor(a-b)
d.ab[is.na(d.ab)] <- 0

# Plot results with SSE annotation
plot.bio3d(abs(d.ab), resno=pdb, sse=pdb, typ="h", xlab="Residue No.", ylab = "Difference Angle")

#DISTANCE MATRIX

a <- dm.xyz(a.xyz)
b <- dm.xyz(b.xyz)

plot.dmat( (a - b), nlevels=10, grid.col="gray", xlab="1tag", ylab="1tnd")