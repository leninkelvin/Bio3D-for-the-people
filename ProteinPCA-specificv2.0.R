library(bio3d)
library(httr)
library(dplyr)
library(openxlsx)
library(DT)


# Download some example PDB files

{ids <- c('1a28',
          '1e3k',
          '1sqn',
          '1sr7',
          '1zuc',
          '2c7a',
          '2ovh',
          '2ovm',
          '2w8y',
          '3d90',
          '3g8o',
          '3hq5',
          '3kba',
          '3zr7',
          '3zra',
          '3zrb',
          '4a2j',
          '4apu',
          '4oar',
          '5cc0'
)} # All PDB

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
         ) # Selected PDBs and chains

raw.files <- get.pdb(ids)

anno <- pdb.annotate(ids)
datatable(anno)

write.xlsx(anno, 'anno.xlsx') # Annotations are manually introduced. Some are resolution, coregulator and ligand type
annotation <- cbind(anno, color)
datatable(annotation)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)

## Calculate RMSD
rmsd(pdbs, fit=TRUE)

core <- core.find(pdbs)

# See Figure 3.
col=rep("black", length(core$volume))
col[core$volume<2]="pink"; col[core$volume<1]="red"
plot(core, col=col)

core.inds <- print(core, vol=1.0)

write.pdb(xyz=pdbs$xyz[1,core.inds$xyz], file="quick_core.pdb")

xyz <- pdbfit( pdbs, core.inds )

rd <- rmsd(xyz)
heatmap(rd, scale = "column")
write.table(rd, file = "RMSD_matrix.tab", sep = "\t", row.names = T)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")

		# Save as PNG
			png("RMSDheatmap.png", width=1200, height=800)
			heatmap(rd, scale = "column")
		dev.off()

# RMSD clustering
hc.rd <- hclust(as.dist(rd))

pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd, k = 4, colors=annotation[, "coregulator"], labels=pdbs$id, cex=1.2,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

	# Save as PNG
		png("RMSD_cluster.png", width=1200, height=800)
		hclustplot(hc.rd, k = 4, colors=annotation[, "coregulator"], labels=pdbs$id, cex=1.2, 
		ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
    dev.off()

# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
pc.xray

plot(pc.xray, col=annotation[, "coregulator"])

    # Save as PNG
    png("PCAx4.png", width=1200, height=800)
    plot(pc.xray, col=annotation[, "coregulator"])
    dev.off()

hc <- hclust(dist(pc.xray$z[,1:2]))
grps <- cutree(hc, h=30)
cols <- c("red", "green", "blue")[grps]
plot(pc.xray, pc.axes=1:2, col=annotation[, "coregulator"])

      # Save as PNG
      png("PCAcluster.png", width=1200, height=800)
      plot(pc.xray, pc.axes=1:2, col=annotation[, "coregulator"])
      dev.off()

# Dendrogram plot
names(cols) <- pdbs$id
hclustplot(hc, col=annotation[, "coregulator"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)

      # Save as PNG
      png("PCA_dendogram.png", width=1200, height=800)
      hclustplot(hc, col=annotation[, "coregulator"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)
      dev.off()




