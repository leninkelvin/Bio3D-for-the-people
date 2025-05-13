library(bio3d)
library(httr)

##WORKING WITH MULTIPLE PDBs

# Download some example PDB files
ids <- c('1a28_A',
         '1e3k_A',
         '1sqn_A',
         '1sr7_A',
         '1zuc_A',
         '2ovh_A',
         '2ovm_A',
         '2w8y_A',
         '3d90_A',
         '3g8o_A',
         '3hq5_A',
         '3kba_A',
         '3zr7_A',
         '3zra_A',
         '3zrb_A',
         '4a2j_A',
         '4apu_A',
         '4oar_A'
)

# Extract and align sequences
#pdb <- read.pdb("myfile.pdb")
#pdb <- read.pdb("/path/to/my/data/myfile.pdb")

raw.files <- get.pdb(ids)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

## Calculate RMSD
rmsd(pdbs, fit=TRUE)

# Extract and align sequences
pdbs <- pdbaln(files)

core <- core.find(pdbs)

# See Figure 3.
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)

core.inds <- print(core, vol=1.0) # good place to check if the alignment worked

write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core.pdb")

xyz <- pdbfit( pdbs, core.inds, outpath="corefit_structures" )

      ##-- Compare to fitting on all equivalent positions
      xyz2 <- pdbfit(pdbs)
      
      ## Note that overall RMSD will be higher but RMSF will
      ##  be lower in core regions, which may equate to a
      ##  'better fit' for certain applications
      gaps <- gap.inspect(pdbs$xyz)
      rmsd(xyz[,gaps$f.inds])
      rmsd(xyz2[,gaps$f.inds])
      
      plot(rmsf(xyz[,gaps$f.inds]), typ="l", col="blue", ylim=c(0,9))
      points(rmsf(xyz2[,gaps$f.inds]), typ="l", col="red")

rd <- rmsd(xyz)
heatmap(rd, scale = "column")
write.table(rd, file = "RMSD_matrix.tab", sep = "\t", row.names = T)

# Save as PNG
png("RMSDheatmap.png", width=1200, height=800)
heatmap(rd, scale = "column")
dev.off()

hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd, k=5, labels=pdbs$id, cex=1.2,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

      # Save as PNG
      png("RMSD_cluster.png", width=1200, height=800)
      hclustplot(hc.rd, k=5, labels=pdbs$id, cex=1.2,
                 ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
      dev.off()

#RMSF
# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# Tailor the PDB structure to exclude gap positions for SSE annotation
pdb = read.pdb("1a28")
id = grep("1a28_A", pdbs$id)
ref.pdb = trim.pdb(pdb, inds=atom.select(pdb, resno = pdbs$resno[id, gaps.res$f.inds]))

# Plot RMSF with SSE annotation and labeled with residue numbers (Figure 8.)
rf <- rmsf(xyz[, gaps.pos$f.inds])
plot.bio3d(rf, resno=ref.pdb, sse=ref.pdb, ylab="RMSF (Å)", xlab="Residue No.", typ="l")

    # Save as PNG
    png("RMSF.png", width=1200, height=800)
    plot.bio3d(rf, resno=ref.pdb, sse=ref.pdb, ylab="RMSF (Å)", xlab="Residue No.", typ="l")
    dev.off()

tor <- torsion.pdb(pdb)

# Basic Ramachandran plot (Figure 9)
plot(tor$phi, tor$psi, xlab="phi", ylab="psi")

    # Save as PNG
    png("torsions.png", width=1200, height=800)
    plot(tor$phi, tor$psi, xlab="phi", ylab="psi")
    dev.off()

# Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
pc.xray

plot(pc.xray)
plot(pc.xray, pc.axes=)

    # Save as PNG
    png("PCAx4.png", width=1200, height=800)
    plot(pc.xray, pc.axes=)
    dev.off()

#plot(pc.xray, pc.axes=1:2)

# Left-click on a point to label and right-click to end
#identify(pc.xray$z[,1:2], labels=basename.pdb(pdbs$id))

par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC1")
plot.bio3d(pc.xray$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC2")
plot.bio3d(pc.xray$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC3")

# See Figure 15
mktrj.pca(pc.xray, pc=1, file="pc1.pdb")
mktrj.pca(pc.xray, pc=2, file="pc2.pdb")
mktrj.pca(pc.xray, pc=3, file="pc3.pdb")

hc <- hclust(dist(pc.xray$z[,1:2]))
grps <- cutree(hc, h=4)
cols <- c("red", "green", "blue")[grps]
plot(pc.xray, pc.axes=1:2, col=cols)

      # Save as PNG
      png("PCAcluster.png", width=1200, height=800)
      plot(pc.xray, pc.axes=1:2, col=cols)
      dev.off()

# Dendrogram plot
names(cols) <- pdbs$id
hclustplot(hc, colors=cols, ylab="Distance in PC Space", main="PC1-2", cex=0.9, fillbox=FALSE, labels = pdbs$id)

      # Save as PNG
      png("PCA_dendogram.png", width=1200, height=800)
      hclustplot(hc, colors=cols, ylab="Distance in PC Space", main="PC1-2", cex=0.9, fillbox=FALSE)
      dev.off()
