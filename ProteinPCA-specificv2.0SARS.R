library(bio3d)
library(httr)
library(dplyr)
library(openxlsx)

# Download some example PDB files
ids <- c( "8C3V_R", "8BCZ_R", "8BBO_R", "7WRO_R", "7WRY_R", "7WRJ_R", 
          "7Y0W_R", "7Y0V_R", "7WRZ_R", "8HWT_A", "8HWS_A", "7X1N_M",
          "7E8F_R", "7CHF_R", "7CH4_R", "7EY0_R", "7EYA_R", "7ZF5_A",
          "7Q0H_E", "7PS6_E", "8SDF_A", "8SDG_C", "8SIT_A", "8D8Q_A", "7QNY_E",
          "7DX4_E", "7MMO_C", "7ZFB_E", "7ZR8_E",
          "7ZRN_E", "7TAS_E", "7RAL_B", "7M7W_R",  
          "8ERQ_A", "7JX3_E", "8S9G_E", "8FXC_E", "8FXB_E", "7XSW_R", 
          "8DW3_B", "8DW2_B", "8IX3_G", "7WTG_E",
          "7WTH_E", "7WTJ_E", "7YD1_C", "7WEF_E"
         )
raw.files <- get.pdb(ids)

anno <- pdb.annotate(ids)
View(anno)

write.xlsx(anno, 'anno.xlsx')
annotation <- cbind(anno, color)

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
hclustplot(hc.rd, k = 4, colors=annotation[, "color"], labels=pdbs$id, cex=1.2,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

	# Save as PNG
		png("RMSD_cluster.png", width=1200, height=800)
		hclustplot(hc.rd, k = 4, colors=annotation[, "color"], labels=pdbs$id, cex=1.2, 
		ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
    dev.off()

# Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

# Do PCA
pc.xray <- pca.xyz(xyz[, gaps.pos$f.inds])
pc.xray

plot(pc.xray, col=annotation[, "color"])

    # Save as PNG
    png("PCAx4.png", width=1200, height=800)
    plot(pc.xray, col=annotation[, "color"])
    dev.off()

hc <- hclust(dist(pc.xray$z[,1:2]))
grps <- cutree(hc, h=30)
cols <- c("red", "green", "blue")[grps]
plot(pc.xray, pc.axes=1:2, col=annotation[, "color"])

      # Save as PNG
      png("PCAcluster.png", width=1200, height=800)
      plot(pc.xray, pc.axes=1:2, col=annotation[, "color"])
      dev.off()

# Dendrogram plot
names(cols) <- pdbs$id
hclustplot(hc, col=annotation[, "color"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)

      # Save as PNG
      png("PCA_dendogram.png", width=1200, height=800)
      hclustplot(hc, col=annotation[, "color"], ylab="Distance in PC Space", main="PC1-2", cex=0.9)
      dev.off()


