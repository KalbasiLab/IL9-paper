###################
Pre-processing
##################
library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(tidyverse)
library(Seurat)
library(scCustomize)
library(Azimuth)

#load the data
IL9Tx_IL9R_AQ<- Read10X(data.dir = "~/Desktop/IL9Tx_IL9R_AQ/filtered_feature_bc_matrix") 
IL9Tx_IL9R_PR<- Read10X(data.dir = "~/Desktop/IL9Tx_IL9R_PR/filtered_feature_bc_matrix") 
IL9Tx_IL9R_WT<- Read10X(data.dir = "~/Desktop//IL9Tx_IL9R_WT/filtered_feature_bc_matrix")

#CRATE Seurat OBJECT
srat_IL9Tx_IL9R_AQ<- CreateSeuratObject(IL9Tx_IL9R_AQ,project = "IL9Tx_IL9R_AQ", min.cells = 3, min.features = 200)
srat_IL9Tx_IL9R_PR <- CreateSeuratObject(IL9Tx_IL9R_PR,project = "IL9Tx_IL9R_PR", min.cells = 3, min.features = 200)
srat_IL9Tx_IL9R_WT<-CreateSeuratObject(IL9Tx_IL9R_WT,project = "IL9Tx_IL9R_WT", min.cells = 3, min.features = 200)

#Add metadata 
srat_IL9Tx_IL9R_AQ <- AddMetaData(srat_IL9Tx_IL9R_AQ ,metadata = "IL9Tx_IL9R_AQ", col.name = "Treatment")
srat_IL9Tx_IL9R_PR <- AddMetaData(srat_IL9Tx_IL9R_PR ,metadata = "IL9Tx_IL9R_PR", col.name = "Treatment")
srat_IL9Tx_IL9R_WT <- AddMetaData(srat_IL9Tx_IL9R_WT ,metadata = "IL9Tx_IL9R_WT", col.name = "Treatment")

#Check for the mitochondria genes
grep ("^mt-", rownames(srat_IL9Tx_IL9R_AQ[["RNA"]]),value = T)
grep ("^mt-", rownames(srat_IL9Tx_IL9R_PR[["RNA"]]),value = T)
grep ("^mt-", rownames(srat_IL9Tx_IL9R_WT[["RNA"]]),value = T)

#MT, percent hemoglobin genes, no. of RNA molecules and number of features
srat_IL9Tx_IL9R_AQ[["percent.mt"]]  <- PercentageFeatureSet(srat_IL9Tx_IL9R_AQ, pattern = "^mt-")
srat_IL9Tx_IL9R_PR[["percent.mt"]]  <- PercentageFeatureSet(srat_IL9Tx_IL9R_PR, pattern = "^mt-")
srat_IL9Tx_IL9R_WT[["percent.mt"]]  <- PercentageFeatureSet(srat_IL9Tx_IL9R_WT, pattern = "^mt-")

#Plot
VlnPlot(srat_IL9Tx_IL9R_AQ, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 4)
VlnPlot(srat_IL9Tx_IL9R_PR, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 4)
VlnPlot(srat_IL9Tx_IL9R_WT, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 4)

#filtering
plot_srat_IL9Tx_IL9R_AQ <- FeatureScatter(srat_IL9Tx_IL9R_AQ, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_srat_IL9Tx_IL9R_AQ
srat_IL9Tx_IL9R_AQ <- subset(srat_IL9Tx_IL9R_AQ, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 15)

plot_srat_IL9Tx_IL9R_PR <- FeatureScatter(srat_IL9Tx_IL9R_PR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_srat_IL9Tx_IL9R_PR
srat_IL9Tx_IL9R_PR <- subset(srat_IL9Tx_IL9R_PR, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 10)

plot_srat_IL9Tx_IL9R_WT <- FeatureScatter(srat_IL9Tx_IL9R_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_srat_IL9Tx_IL9R_WT
srat_IL9Tx_IL9R_WT <- subset(srat_IL9Tx_IL9R_WT, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 15)

#Merge these into one object
merged_obj<- merge(x=srat_IL9Tx_IL9R_AQ ,y=list(srat_IL9Tx_IL9R_PR,srat_IL9Tx_IL9R_WT))
merged_obj <- JoinLayers(merged_obj)

#split 
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f=merged_obj$Treatment)

#Perform analysis without integration 
#https://satijalab.org/seurat/articles/integration_introduction.html
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
elbow <- ElbowPlot(merged_obj, ndims = 30)
elbow
merged_obj <- RunHarmony(merged_obj, group.by.vars = "Treatment") #this should be saved as another slot

merged_obj <- FindNeighbors(merged_obj, dims = 1:17, reduction = "pca") #PCA  USE harmony reductions if they are processed separately (to remove bacth effetcs)
#merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "pca") 

#resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6)
#nochemo_pre <- FindClusters(nochemo_pre, resolution = resolution.range)

merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "merged_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:17, reduction = "pca") #pca 
#merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca") #pca 

merged_obj<-RunTSNE(merged_obj, reduction = "pca", dims = 1:17) 
#merged_obj<-RunTSNE(merged_obj, reduction = "pca", dims = 1:30) 

# can change dims.use to see if it affects the clustering
DimPlot(merged_obj,  reduction = "umap", group.by = c("Treatment", "seurat_clusters"))
DimPlot(merged_obj, reduction = "tsne", split.by = c("Treatment"))
DimPlot(merged_obj, reduction = "umap", split.by = c("Treatment"), label=T)

#This is 0.5 resolution object this is the object we used 
#saveRDS(merged_obj, file = "~/Desktop/Kayla_scData/merged_obj.rds") 
merged_obj <- readRDS("~/Desktop/Kayla_scData/merged_obj.rds")
merged_obj <- JoinLayers(merged_obj)
#now subset the any other group to only include these three.
Idents(merged_obj) <-"Treatment"
Idents(merged_obj_subset) <-"seurat_clusters"
DefaultAssay(merged_obj) <- "RNA"
#Subset
merged_obj_subset <- subset(x = merged_obj, idents = c("IL9Tx_IL9R_AQ", "IL9Tx_IL9R_PR","IL9Tx_IL9R_WT"))
table(Idents(merged_obj_subset))

DefaultAssay(merged_obj) <- "RNA"
FeaturePlot(merged_obj, features = "Ptprc", split.by = "Treatment")
FeaturePlot(merged_obj, features = "Thy1")

DimPlot_scCustom(seurat_object = merged_obj,colors_use = DiscretePalette_scCustomize(num_colors = 24, 
                                                                                     palette = "ditto_seq"))                                                                                                                                                                                 
DimPlot_scCustom(seurat_object = merged_obj,split.by = "Treatment", colors_use = DiscretePalette_scCustomize(num_colors = 24, 
                                                                                                             palette = "ditto_seq"))                                                                                                                                                                                 



#plotting

#Changed the colors 
cluster_colors_treatment <- c("Tem/rm" = "#E64B35FF",  
                              "T early activated" = "#FF7F00",
                              "Teff-1/prolif" = "#984EA3",
                              "Teff-2/prolif" = "#4DBBD5FF",
                              "Teff-3/pex/naïve(mix)"="#42B540FF","Teff-4/naïve/em(mix)"="#CCAA7A",
                              "Teff-5(Stat1/Gzmb/Prf1)"="#7E6148FF","Tscm-like(Tcf7)"="#7D3434","Tmyeloid-like"="#79AF97FF","Tprolif"="#3C5488FF")

##Changed the colors --MAKE these in the barplot and other ridgeplots
#Plotting

p1 <- DimPlot(merged_obj_subset, reduction = "umap", split.by = "Treatment", label = FALSE, cols = cluster_colors_treatment,
              pt.size = 0.5)+
  #ggtitle("Clusters") +
  scale_color_manual(values = cluster_colors_treatment) + 
  theme(axis.title.x = element_text(color="black", size=28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=28, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(color="black", size=28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y = element_text(color="black", size=28, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=28, color="black"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) # Increase overall plot margins

p1

ggsave("~/Desktop/Figure1-IL9/SPLIT_UMAP.pdf", width = 36, height = 18, units = "cm") # Slightly increased dimensions


##Changed the colors --MAKE these in the barplot and other ridgeplots
#Plotting

p2 <- DimPlot(merged_obj_subset, reduction = "umap", label = FALSE, cols = cluster_colors_treatment,pt.size = 0.5)+
  #ggtitle("Clusters") +
  scale_color_manual(values = cluster_colors_treatment) + 
  theme(axis.title.x = element_text(color="black", size=28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=28, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(color="black", size=28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y = element_text(color="black", size=28, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=28, color="black"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) # Increase overall plot margins

p2
ggsave("~/Desktop/Figure1-IL9/UMAP.pdf", width = 28, height = 18, units = "cm") # Slightly increased dimensions

########
#Color the UMAP with the three colors

#pstat3
cluster_colors_treatment <- c("IL9Tx_IL9R_WT" = "blue",
                              "IL9Tx_IL9R_PR" = "#F39B7FFF",
                              "IL9Tx_IL9R_AQ" = "#91D1C2FF")
Idents(merged_obj_subset) <-"Treatment"

p3 <- DimPlot(merged_obj_subset,reduction = "umap", label = FALSE, cols = cluster_colors_treatment,pt.size = 0.75)+
  #ggtitle("Clusters") +
  scale_color_manual(values = cluster_colors_treatment) + 
  theme(axis.title.x = element_text(color="black", size=28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=28, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(color="black", size=28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.y = element_text(color="black", size=28, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=28, color="black"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) # Increase overall plot margins

p3
ggsave("~/Desktop/Figure1-IL9/UMAP_bytreatment.pdf", width = 36, height = 18, units = "cm") # Slightly increased dimensions



#####Barplot
library(plot1cell)
library(plot1cell)
library(ggplot2)
library(speckle)
library(limma)
library(ggplot2)

plotCellTypeProps <- function(x = NULL, clusters = NULL, sample = NULL, custom_colors = NULL)
{
  if (is.null(x) & is.null(sample) & is.null(clusters))
    stop("Please provide either a SingleCellExperiment object or Seurat
          object with required annotation metadata, or explicitly provide
          clusters and sample information")
  
  if ((is.null(clusters) | is.null(sample)) & !is.null(x)) {
    # Extract cluster, sample, and group info from SCE object
    if (is(x, "SingleCellExperiment"))
      y <- .extractSCE(x)
    
    # Extract cluster, sample, and group info from Seurat object
    if (is(x, "Seurat"))
      y <- .extractSeurat(x)
    
    clusters <- y$clusters
    sample <- y$sample
  }
  
  prop.list <- getTransformedProps(clusters, sample)
  
  Proportions <- as.vector(t(prop.list$Proportions))
  Samples <- rep(colnames(prop.list$Proportions), nrow(prop.list$Proportions))
  Clusters <- rep(rownames(prop.list$Proportions),
                  each = ncol(prop.list$Proportions))
  
  plotdf <- data.frame(Samples = Samples, Clusters = Clusters,
                       Proportions = Proportions)
  # Relevel the Samples factor
  plotdf$Samples <- factor(plotdf$Samples, levels = c("IL9Tx_IL9R_WT", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_AQ"))
  
  gg <- ggplot(plotdf, aes(x = Samples, y = Proportions, fill = Clusters)) +
    geom_bar(stat = "identity", color = "black", width = 0.5) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20,colour = "black"),
          axis.text.y = element_text(size = 24, colour = "black"),
          axis.title = element_text(size = 16,colour = "black"),
          legend.text = element_text(size = 16,colour = "black"),
          legend.title = element_text(size = 16,colour = "black"),
          plot.title = element_text(size = 16, hjust = 0,colour = "black")) +
    labs(fill = "Clusters")
  
  return(gg)
}

# Your custom colors
#Changed the colors 
cluster_colors_treatment <- c("Tem/rm" = "#E64B35FF",  
                              "T early activated" = "#FF7F00",
                              "Teff-1/prolif" = "#984EA3",
                              "Teff-2/prolif" = "#4DBBD5FF",
                              "Teff-3/pex/naïve(mix)"="#42B540FF","Teff-4/naïve/em(mix)"="#CCAA7A",
                              "Teff-5(Stat1/Gzmb/Prf1)"="#7E6148FF","Tscm-like(Tcf7)"="#7D3434","Tmyeloid-like"="#79AF97FF","Tprolif"="#3C5488FF")

#Create the plot
p <- plotCellTypeProps(clusters = Idents(merged_obj_subset), 
                       sample = merged_obj_subset$Treatment, 
                       custom_colors = cluster_colors_treatment)
# Display the plot
print(p)

ggsave("~/Desktop/Figure1-IL9/barplot_byTreatment.pdf", width = 20, height = 17, units = "cm")


#Collapsed version
merged_obj_subset <- RenameIdents(merged_obj_subset, `0` = "Tem/rm", `1` = "Teff-1/activated", `2` = "Teff-2/prolif",
                                  `3` = "Teff-3/prolif", `4` = "Teff-4/pex/naïve(mix)", `5` = "Teff-5/naïve/em(mix)", `6` = "Teff-6(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")
#rename again
merged_obj_subset_collaped <- RenameIdents(merged_obj_subset, `Tem/rm` = "Tem/rm", `Teff-1/activated` = "Teff", `Teff-2/prolif` = "Teff",
                                           `Teff-3/prolif` = "Teff", `Teff-4/pex/naïve(mix)` = "Teff", `Teff-5/naïve/em(mix)` = "Teff", `Teff-6(Stat1/Gzmb/Prf1)` = "Teff", `Tscm-like(Tcf7)` = "Tscm-like(Tcf7)",`Tmyeloid-like`="Tmyeloid-like",`Tprolif`="Tprolif")

# Step 1: Rename identities using character strings
merged_obj_subset_collaped <- RenameIdents(merged_obj_subset, 
                                           `Tem/rm` = "Tem/rm", 
                                           `Teff-1/activated` = "Teff", 
                                           `Teff-2/prolif` = "Teff",
                                           `Teff-3/prolif` = "Teff", 
                                           `Teff-4/pex/naïve(mix)` = "Teff", 
                                           `Teff-5/naïve/em(mix)` = "Teff", 
                                           `Teff-6(Stat1/Gzmb/Prf1)` = "Teff", 
                                           `Tscm-like(Tcf7)` = "Tscm-like(Tcf7)",
                                           `Tmyeloid-like` = "Tmyeloid-like",
                                           `Tprolif` = "Tprolif")

# Step 2: Create a named vector of expressions for plotting
label_expressions <- c(
  Teff = expression(T[eff]),
  `Tem/rm` = expression(T[em/rm]),
  `Tscm-like(Tcf7)` = expression(T[scm-like(Tcf7)]),
  `Tmyeloid-like` = expression(T[myeloid-like]),
  Tprolif = expression(T[prolif])
)


library(ggplot2)
ggplot(data, aes(x = Idents(merged_obj_subset_collaped), y = some_value)) +
  geom_boxplot() +
  scale_x_discrete(labels = label_expressions) +
  theme_minimal()

###merged_obj_subset_collaped
plotCellTypeProps(clusters = Idents(merged_obj_subset_collaped), sample = merged_obj_subset_collaped$Treatment) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = label_expressions)
ggtitle("cell type proportions") +
  theme(plot.title = element_text(size = 18, hjust = 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))

