library(bio3d)
library(httr)

##WORKING WITH TRAJECTORIES

#dcdfile <- "C:/Users/acher/Downloads/examples/hivp.dcd"
#pdbfile <- "C:/Users/acher/Downloads/examples/hivp.pdb"
#dcd <- read.dcd(dcdfile)
#pdb <- read.pdb(pdbfile)

dcdfile <- "C:/Users/acher/Downloads/examples/2LAO_K.dcd" #Windows
pdbfile <- "C:/Users/acher/Downloads/examples/2LAO_K.pdb" #Windows

#dcdfile <- "/Users/kelvin/Downloads/2LAO_K.dcd" #unix
#pdbfile <- "/Users/kelvin/Downloads/2LAO_K.pdb" #unix

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

sse <- dssp(pdb) 

helix_inds <- which(sse$sse == "H")

ca.inds <- atom.select(helix_inds, elety="CA")
#ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)

rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")

pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )

hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)

plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")

p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2.pdb")

write.ncdf(p1, "trj_pc1.nc")

cij<-dccm(xyz[,ca.inds$xyz])
plot(cij)

# View the correlations in pymol
pymol.dccm(cij, pdb, type="launch")

