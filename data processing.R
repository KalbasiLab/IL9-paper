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







####Trajectory Analysis
#### Create a Monocle CDS Object
#this is the 0.5 resolution object
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(Seurat)
set.seed(1)

merged_obj <- readRDS("~/Desktop/Kayla_scData/merged_obj.rds")
merged_obj <- JoinLayers(merged_obj)
DefaultAssay(merged_obj) <- "RNA"

DimPlot(merged_obj, reduction = "umap", label = FALSE)

Idents(merged_obj) <-"Treatment"
DefaultAssay(merged_obj) <- "RNA"
merged_obj_subset <- subset(x = merged_obj, idents = c("IL9Tx_IL9R_AQ", "IL9Tx_IL9R_PR","IL9Tx_IL9R_WT"))

Idents(merged_obj_subset) <-"seurat_clusters"
#Monocle3
# Project PC dimensions to whole data set
gene_annotation <- as.data.frame(rownames(merged_obj_subset@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(merged_obj_subset@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(merged_obj_subset[["RNA"]]$counts@Dimnames[[2]],
                               row.names = merged_obj_subset[["RNA"]]$counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
New_matrix <- merged_obj_subset[["RNA"]]$counts
New_matrix <- New_matrix[rownames(merged_obj_subset@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

####
cds_from_seurat <- new_cell_data_set(expression_matrix,cell_metadata = cell_metadata,gene_metadata = gene_annotation) #created the monocle3 object

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- merged_obj_subset@active.ident
names(list_cluster) <- merged_obj_subset[["RNA"]]$counts@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- merged_obj_subset@reductions[["umap"]]@cell.embeddings

#Then perform pseudotime analysis:
cds_from_seurat <- learn_graph(cds_from_seurat)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=4, cell_size=0.5,group_label_size=5)

#calculate the psudotime for each cell after defining the root cells 
#cds_from_seurat <-  order_cells(cds_from_seurat, reduction_method = "UMAP", root_cells = colnames(cds_from_seurat[, clusters(cds_from_seurat) == "T-stem-like (Tcf7+)"]))
#ps_tim <- pseudotime(cds_from_seurat)

cds_from_seurat <-  order_cells(cds_from_seurat,reduction_method = "UMAP")

plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, #Branch points
           group_label_size = 1,
           graph_label_size=4,label_roots = TRUE,cell_size = 0.8)+
  theme(
    axis.text = element_text(size = 19, colour = "black"),  # Increase axis text size
    axis.title = element_text(size = 17,colour = "black"),  # Increase axis title size
    legend.text = element_text(size = 17,colour = "black"),  # Increase legend text size
    legend.title = element_text(size = 17,colour = "black")  # Increase legend title size
  )

#ggsave("pseudotime_umap_new.pdf", width = 13, height = 12, units = "cm")


#scale_colour_gradient2(low = "red", mid = "white",high = "blue", midpoint = 3, space = "Lab",na.value = "#F5F5F500", guide = "colourbar", aesthetics = "colour")

#ggsave("pseudotime_umap.pdf", width = 18, height = 12, units = "cm")
#ggsave("~/UMAP_monocle3.pdf", width = 18, height = 12, units = "cm")

#Extract the pseudotime value and add it to the seurat
merged_obj_subset <- AddMetaData(
  object = merged_obj_subset,
  metadata = cds_from_seurat@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime_scores"
)

#Plot this as a violin plot
FeaturePlot(merged_obj_subset, "pseudotime_scores", pt.size = 0.1) & scale_color_viridis_c()

VlnPlot(merged_obj, features = "pseudotime_scores")
Idents(merged_obj) <- merged_obj$Treatment
VlnPlot(merged_obj, features = c("pseudotime_scores"), pt.size = 0, ncol = 1,add.noise = FALSE)

c("Tem/rm" = "#E64B35FF",  
  "T early activated" = "#FF7F00",
  "Teff-1/prolif" = "#984EA3",
  "Teff-2/prolif" = "#4DBBD5FF",
  "Teff-3/pex/naïve(mix)"="#42B540FF","Teff-4/naïve/em(mix)"="#CCAA7A",
  "Teff-5(Stat1/Gzmb/Prf1)"="#7E6148FF","Tscm-like(Tcf7)"="#7D3434","Tmyeloid-like"="#79AF97FF","Tprolif"="#3C5488FF")

Idents(merged_obj_subset) <-"seurat_clusters"

merged_obj_subset <- RenameIdents(merged_obj_subset, `0` = "Tem/rm", `1` = "T early activated", `2` = "Teff-1/prolif",
                                  `3` = "Teff-2/prolif", `4` = "Teff-3/pex/naïve(mix)", `5` = "Teff-4/naïve/em(mix)", `6` = "Teff-5(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")
#pLOTTING
RidgePlot(merged_obj_subset, features = "pseudotime_scores") +  
  scale_fill_manual(values = c("Tem/rm" = "#E64B35FF",  
                               "T early activated" = "#FF7F00",
                               "Teff-1/prolif" = "#984EA3",
                               "Teff-2/prolif" = "#4DBBD5FF",
                               "Teff-3/pex/naïve(mix)"="#42B540FF","Teff-4/naïve/em(mix)"="#CCAA7A",
                               "Teff-5(Stat1/Gzmb/Prf1)"="#7E6148FF","Tscm-like(Tcf7)"="#7D3434","Tmyeloid-like"="#79AF97FF","Tprolif"="#3C5488FF"))+
  theme_minimal() +
  labs(title = "Pseudotime Scores by Cluster",
       y = "Cluster Annotation",
       x = "Pseudotime Scores") +
  theme(axis.text.x = element_text(hjust = 1, color="black", size=17),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(color="black", size=17),
        axis.title.x = element_text(color="black", size=17),
        axis.text.y = element_text(color="black", size=17),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 17, hjust = 0.5))

ggsave("pseudotime_RidgePlot_by_clusters.pdf", width = 24, height = 12, units = "cm")



############
Idents(merged_obj_subset) <- "Treatment"
selected_conditions <- c("IL9Tx_IL9R_AQ", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_WT")
# Ensure the condition is a factor
merged_obj_subset$Treatment <- factor(merged_obj_subset$Treatment, levels = selected_conditions)
# Set the Idents of the Seurat object
Idents(merged_obj_subset) <- "Treatment"


RidgePlot(merged_obj_subset, features = "pseudotime_scores")+
  scale_fill_manual(values = c("#91D1C2FF", "#F39B7FFF", "blue")) +
  theme_minimal() +
  labs(title = "Pseudotime Scores by Treatment",
       y = "Treatment",
       x = "Pseudotime Scores") +
  theme(axis.text.x = element_text(hjust = 1,color="black", size=14),
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.title.y = element_text(color="black", size=14),
        axis.title.x = element_text(color="black", size=14),
        axis.text.y= element_text(color="black", size=14),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.background = element_blank(),   # Remove plot background
        panel.border = element_rect(color = "black", fill = NA,size = 1),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 14, hjust = 0.5)) 



#reorder:
# Define the order for plotting
plot_order <- c("IL9Tx_IL9R_AQ", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_WT")

# Define the order for the legend
legend_order <- c("IL9Tx_IL9R_WT", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_AQ")

# Ensure the condition is a factor with the plot order
merged_obj_subset$Treatment <- factor(merged_obj_subset$Treatment, levels = plot_order)

# Set the Idents of the Seurat object
Idents(merged_obj_subset) <- "Treatment"

# Create the RidgePlot
RidgePlot(merged_obj_subset, features = "pseudotime_scores") +
  scale_fill_manual(values = c("IL9Tx_IL9R_AQ" = "#91D1C2FF", 
                               "IL9Tx_IL9R_PR" = "#F39B7FFF", 
                               "IL9Tx_IL9R_WT" = "blue"),
                    breaks = legend_order) +
  theme_minimal() +
  labs(title = "Pseudotime Scores by Treatment",
       y = "Treatment",
       x = "Pseudotime Scores") +
  theme(axis.text.x = element_text(hjust = 1, color="black", size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.y = element_text(color="black", size=16),
        axis.title.x = element_text(color="black", size=16),
        axis.text.y = element_text(color="black", size=16),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))


ggsave("pseudotime_RidgePlot_by_treatment.pdf", width = 18, height = 10, units = "cm")


#BOXPLOT
library(Seurat)
library(ggplot2)

DefaultAssay(merged_obj_subset) <- "RNA"
Idents(merged_obj_subset) <- "seurat_clusters"
#Rename the idents
merged_obj_subset <- RenameIdents(merged_obj_subset, `0` = "Tem/rm", `1` = "Teff-1/activated", `2` = "Teff-2/prolif",
                                  `3` = "Teff-3/prolif", `4` = "Teff-4/pex/naïve(mix)", `5` = "Teff-5/naïve/em(mix)", `6` = "Teff-6(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")



merged_obj_subset$cluster_annotation <- Idents(merged_obj_subset)

selected_conditions <- c("IL9Tx_IL9R_WT", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_AQ")

#Extract data from Seurat object for ggplot
data <- FetchData(merged_obj_subset, vars = c("pseudotime_scores", "Treatment","cluster_annotation"))
data_filtered <- data[data$Treatment %in% selected_conditions, ]
# Ensure the condition is a factor
data_filtered$Treatment <- factor(data_filtered$Treatment, levels = selected_conditions)
#Create a boxplot using ggplot2

boxplot <- ggplot(data_filtered, aes(x = Treatment, y = pseudotime_scores, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("blue", "#F39B7FFF", "#91D1C2FF")) +
  theme_minimal() +
  labs(title = "Pseudotime Scores by Treatment",
       x = "Treatment",
       y = "Pseudotime Scores") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black", size=10),
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.text.y= element_text(color="black", size=10),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.background = element_blank(),   # Remove plot background
        panel.border = element_rect(color = "black", fill = NA))  # Add border around plot area

print(boxplot)
ggsave("pseudotime_boxplot_by_treatment.pdf", width = 12, height = 8, units = "cm")


#By cluster annotation
# Create a boxplot using ggplot2
boxplot <- ggplot(data_filtered, aes(x = cluster_annotation, y = pseudotime_scores, fill = cluster_annotation)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666", "#AD7700",
                               "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF",
                               "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C")) +
  theme_minimal() +
  labs(title = "Pseudotime Scores by cluster annotation",
       x = "Cluster Annotation",
       y = "Pseudotime Scores") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black", size=10),
        panel.background = element_blank(),  # Remove panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        axis.text.y= element_text(color="black", size=10),
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.background = element_blank(),   # Remove plot background
        panel.border = element_rect(color = "black", fill = NA))  # Add border around plot area

print(boxplot)

ggsave("~/Desktop/Figure1-IL9/pseudotime_boxplot_by_cluster-annotation.pdf", width = 12, height = 10, units = "cm")


##pstat1 module vs pseudotime
# Extract module scores
pseudotime_scores <- merged_obj_subset[["pseudotime_scores"]]
pstat1_scores <- score_pstat1[["pstat11"]]

#Combine scores into a data frame
score_df <- data.frame(
  pseudotime = pseudotime_scores,
  pSTAT1 = pstat1_scores
)

#Calculate correlation
correlation <- cor(score_df$pseudotime_scores, score_df$pstat11, method = "spearman")

# Print the correlation
print(paste("Correlation between pseudotime and pSTAT1 scores:", correlation))

#Visualize the correlation with a scatter plot
library(ggplot2)
ggplot(score_df, aes(x = pseudotime_scores, y = pstat11)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Correlation between pseudotime and pSTAT1 scores",
       x = "pseudotime score",
       y = "pSTAT1 score") +
  theme_minimal()



###
#TCF7 
VlnPlot(merged_obj_subset, features = c("Tcf7"), pt.size = 0, ncol = 1, add.noise = FALSE, cols = cluster_colors_treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 18),
    axis.text.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", size = 18),
    legend.text = element_text(size = 18, color = "black")
  ) +
  ggtitle("tcf7")

ggsave("Tcf7_by_treatment.pdf", width = 18, height = 12, units = "cm")


#pstat3
Idents(merged_obj_subset) <-"Treatment"

cluster_colors_treatment <- c("IL9Tx_IL9R_WT" = "blue",
                              "IL9Tx_IL9R_PR" = "#F39B7FFF",
                              "IL9Tx_IL9R_AQ" = "#91D1C2FF")

VlnPlot(merged_obj_subset, features = c("Il7r"), pt.size = 0, ncol = 1, add.noise = FALSE, cols = cluster_colors_treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 18),
    axis.text.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", size = 18),
    legend.text = element_text(size = 18, color = "black")
  ) +
  ggtitle("Il7r")

ggsave("Il7r_by_treatment.pdf", width = 18, height = 12, units = "cm")


#IL7R
#by treatment
Idents(merged_obj_subset) <-"seurat_clusters"

merged_obj_subset <- RenameIdents(merged_obj_subset, `0` = "Tem/rm", `1` = "T early activated", `2` = "Teff-1/prolif",
                                  `3` = "Teff-2/prolif", `4` = "Teff-3/pex/naïve(mix)", `5` = "Teff-4/naïve/em(mix)", `6` = "Teff-5(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")


cluster_colors_treatment <- c("Tem/rm" = "#E64B35FF",  
                              "T early activated" = "#FF7F00",
                              "Teff-1/prolif" = "#984EA3",
                              "Teff-2/prolif" = "#4DBBD5FF",
                              "Teff-3/pex/naïve(mix)"="#42B540FF","Teff-4/naïve/em(mix)"="#CCAA7A",
                              "Teff-5(Stat1/Gzmb/Prf1)"="#7E6148FF","Tscm-like(Tcf7)"="#7D3434","Tmyeloid-like"="#79AF97FF","Tprolif"="#3C5488FF")

cluster_colors_treatment <- c("IL9Tx_IL9R_WT" = "blue",
                              "IL9Tx_IL9R_PR" = "#F39B7FFF",
                              "IL9Tx_IL9R_AQ" = "#91D1C2FF")

DefaultAssay(merged_obj_subset) <- "RNA"

VlnPlot(merged_obj_subset, features = c("Tcf7"), pt.size = 0, ncol = 1, add.noise = FALSE, cols = cluster_colors_treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 13),
    axis.text.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    legend.text = element_text(size = 13, color = "black")
  ) +
  ggtitle("Tcf7")

ggsave("Tcf7_by_cluster.pdf", width = 16, height = 10, units = "cm")



#reorder the idents
new_order <- c("IL9Tx_IL9R_WT", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_AQ")
# Reorder the Idents
Idents(merged_obj_subset) <- factor(Idents(merged_obj_subset), levels = new_order)
# Check the new order
table(Idents(merged_obj_subset))


VlnPlot(merged_obj_subset, features = c("Prf1"), pt.size = 0, ncol = 1, add.noise = FALSE, cols = cluster_colors_treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 18),
    axis.text.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", size = 18),
    legend.text = element_text(size = 18, color = "black")
  ) +
  ggtitle("Prf1")

ggsave("prf1_by_treatment.pdf", width = 18, height = 12, units = "cm")



VlnPlot(merged_obj_subset, features = c("Gzmb"), pt.size = 0, ncol = 1, add.noise = FALSE, cols = cluster_colors_treatment) +
  theme(
    axis.title.x = element_text(color = "black", size = 18),
    axis.text.x = element_text(color = "black", size = 18),
    axis.title.y = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", size = 18),
    legend.text = element_text(size = 18, color = "black")
  ) +
  ggtitle("Gzmb")

ggsave("GZMB_by_treatment.pdf", width = 18, height = 12, units = "cm")



