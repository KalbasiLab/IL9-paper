#Date: October 17th
#SCENIC MORE refined analysis

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(readr)

#Load the seurat object
merged_obj <- readRDS("~/Desktop/Kayla_scData/merged_obj.rds")
merged_obj <- JoinLayers(merged_obj)
DefaultAssay(merged_obj) <-"RNA"

#by treatment
Idents(merged_obj) <-"seurat_clusters"

merged_obj <- RenameIdents(merged_obj, `0` = "Tem/rm", `1` = "T early activated", `2` = "Teff-1/prolif",
                                  `3` = "Teff-2/prolif", `4` = "Teff-3/pex/naïve(mix)", `5` = "Teff-4/naïve/em(mix)", `6` = "Teff-5(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")

#new regulons script
regulons <- read_csv("~/Desktop/Kayla_scData/scenic/regulons_new.csv")
head(regulons)

#Load the AUC Matrix
regulonAUC <- importAUCfromText("~/Desktop/Kayla_scData/scenic/auc_mtx_new.csv")  ##New 
Idents(merged_obj) <- "Treatment"
cellInfo <- data.frame(Treatment=Idents(merged_obj))
cellInfo <- as.data.frame(cellInfo)

cellsPerCluster <- split(rownames(cellInfo), cellInfo$Treatment) 

#cellsPerCluster <- as.data.frame(cellsPerCluster)
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = F, scale=T))

temp <- as.data.frame(regulonActivity_byCellType)
temp$tf <- rownames(temp) #This is the relative exrichment score This is achieved using AUCell, which provides an AUC value indicating the relative expression of the genes in the regulon (as a whole) in each cell.

#####
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "Treatment", "RelativeActivity")
topRegulators$Treatment <- factor(as.character(topRegulators$Treatment))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators) #These are the top regulators by treatment group 
#write.csv(topRegulators, "/Users/apal6/Desktop/Kayla_scData/topRegulators_by_treatment_noPruning.csv", row.names = FALSE)



###Integrate scenic analysis in seurat for downstream analysis and plotting
regulonAUC #This is the regulon activity dataframe that we would want to store in the seurat object 
# First, ensure that the cell names in regulonAUC match those in merged_obj
all(colnames(regulonAUC) == colnames(merged_obj))
length(colnames(merged_obj))
#6746 --cell numbers match
# If they don't match, you might need to reorder regulonAUC
regulonAUC <- regulonAUC[, colnames(merged_obj)] 
dim(regulonAUC)
#1109 6746 --cell numbers match
# Convert regulonAUC to a matrix if it's not already
# Try to extract the AUC matrix using the getAUC function
regulonAUC_matrix <- getAUC(regulonAUC)
# Check the dimensions of the resulting matrix
print(dim(regulonAUC_matrix))
#Create a new assay
regulon_assay <- CreateAssayObject(counts = regulonAUC_matrix)

#Add the new assay to the Seurat object
merged_obj[["regulonAUC"]] <- regulon_assay

# Verify that the new assay has been added
print(Assays(merged_obj))

#Drop the extra ident
DefaultAssay(merged_obj) <-"RNA"
Idents(merged_obj) <- "Treatment"
merged_obj_subset <- subset(x = merged_obj, idents = c("IL9Tx_IL9R_AQ", "IL9Tx_IL9R_PR","IL9Tx_IL9R_WT"))
Idents(merged_obj_subset) <- "Treatment"
DefaultAssay(merged_obj_subset) <- "regulonAUC"
#reorder the idents
new_order <- c("IL9Tx_IL9R_WT", "IL9Tx_IL9R_PR", "IL9Tx_IL9R_AQ")
# Reorder the Idents
Idents(merged_obj_subset) <- factor(Idents(merged_obj_subset), levels = new_order)
# Check the new order
table(Idents(merged_obj_subset))

############
#get the ridge plot for "Stat1(+)", "Stat2(+)","Stat3(+)","Stat4(+)","Stat5a(+)","Stat5b(+)"
###########
RidgePlot(object = merged_obj, features = c("Stat1(+)","Stat2(+)","Stat3(+)","Stat4(+)","Stat5a(+)","Stat5b(+)"), assay="regulonAUC")  

#Now perform the analysis on subset object 
cluster_colors_treatment <- c("IL9Tx_IL9R_WT" = "blue",
                              "IL9Tx_IL9R_PR" = "#F39B7FFF",
                              "IL9Tx_IL9R_AQ" = "#91D1C2FF")

RidgePlot(object = merged_obj_subset, features = c("Stat1(+)"), assay="regulonAUC", cols=cluster_colors_treatment)+
  theme(
    axis.title.x = element_text(color = "black", size = 22),
    axis.text.x = element_text(color = "black", size = 22),
    axis.title.y = element_text(color = "black", size = 22),
    axis.text.y = element_text(color = "black", size = 22),
    legend.text = element_text(size = 22, color = "black")
  )

ggsave("~/Desktop/Figure1-IL9/pyscenic/stats1_ridge.pdf", width = 23, height = 14, units = "cm")

#Stat3(+)
RidgePlot(object = merged_obj_subset, features = c("Stat3(+)"), assay="regulonAUC", cols=cluster_colors_treatment)+
  theme(
    axis.title.x = element_text(color = "black", size = 22),
    axis.text.x = element_text(color = "black", size = 22),
    axis.title.y = element_text(color = "black", size = 22),
    axis.text.y = element_text(color = "black", size = 22),
    legend.text = element_text(size = 22, color = "black")
  )
ggsave("~/Desktop/Figure1-IL9/pyscenic/stats3_ridge.pdf", width = 24, height = 14, units = "cm")

#"Stat4(+)"
RidgePlot(object = merged_obj_subset, features = c("Stat4(+)"), assay="regulonAUC", cols=cluster_colors_treatment)+
  theme(
    axis.title.x = element_text(color = "black", size = 22),
    axis.text.x = element_text(color = "black", size = 22),
    axis.title.y = element_text(color = "black", size = 22),
    axis.text.y = element_text(color = "black", size = 22),
    legend.text = element_text(size = 22, color = "black")
  )
ggsave("~/Desktop/Figure1-IL9/pyscenic/stats4_ridge.pdf", width = 23, height = 14, units = "cm")


#"Stat5a(+)"
RidgePlot(object = merged_obj_subset, features = c("Stat5a(+)"), assay="regulonAUC", cols=cluster_colors_treatment)+
  theme(
    axis.title.x = element_text(color = "black", size = 22),
    axis.text.x = element_text(color = "black", size = 22),
    axis.title.y = element_text(color = "black", size = 22),
    axis.text.y = element_text(color = "black", size = 22),
    legend.text = element_text(size = 22, color = "black")
  )
ggsave("~/Desktop/Figure1-IL9/pyscenic/stats5A_ridge.pdf", width = 23, height = 14, units = "cm")




#Put the statistics on top
library(dunn.test)
# For STAT1
stat1_values <- GetAssayData(merged_obj_subset, assay = "regulonAUC", slot = "data")["Stat1(+)",]
kruskal_result_stat1 <- kruskal.test(stat1_values ~ Idents(merged_obj_subset))
dunn_result_stat1 <- dunn.test(stat1_values, g = Idents(merged_obj_subset), method = "holm")
dunn_result_stat1$P
[1] 4.436903e-115  3.762776e-22  5.659261e-65
> dunn_result_stat1$comparisons
[1] "IL9Tx_IL9R_AQ - IL9Tx_IL9R_PR" "IL9Tx_IL9R_AQ - IL9Tx_IL9R_WT" "IL9Tx_IL9R_PR - IL9Tx_IL9R_WT"

#pstat2

stat2_values <- GetAssayData(merged_obj_subset, assay = "regulonAUC", slot = "data")["Stat2(+)",]
kruskal_result_stat2 <- kruskal.test(stat2_values ~ Idents(merged_obj_subset))
dunn_result_stat2 <- dunn.test(stat2_values, g = Idents(merged_obj_subset), method = "holm")
dunn_result_stat2$P
[1] 9.429813e-70 2.292133e-10 2.616608e-97
dunn_result_stat2$comparisons
"IL9Tx_IL9R_AQ - IL9Tx_IL9R_PR" "IL9Tx_IL9R_AQ - IL9Tx_IL9R_WT" "IL9Tx_IL9R_PR - IL9Tx_IL9R_WT"


#pstat3
stat3_values <- GetAssayData(merged_obj_subset, assay = "regulonAUC", slot = "data")["Stat3(+)",]
kruskal_result_stat3 <- kruskal.test(stat3_values ~ Idents(merged_obj_subset))
dunn_result_stat3 <- dunn.test(stat3_values, g = Idents(merged_obj_subset), method = "holm")
#get the stat values
dunn_result_stat3$P
1.112587e-89 4.441314e-78 2.987546e-20
dunn_result_stat3$comparisons
"IL9Tx_IL9R_AQ - IL9Tx_IL9R_PR" "IL9Tx_IL9R_AQ - IL9Tx_IL9R_WT" "IL9Tx_IL9R_PR - IL9Tx_IL9R_WT"


#pstat4
stat4_values <- GetAssayData(merged_obj_subset, assay = "regulonAUC", slot = "data")["Stat4(+)",]
kruskal_result_stat4 <- kruskal.test(stat4_values ~ Idents(merged_obj_subset))
dunn_result_stat4 <- dunn.test(stat4_values, g = Idents(merged_obj_subset), method = "holm")
#get the stat values
dunn_result_stat4$P
1.222297e-292 2.806724e-272  7.937750e-58

dunn_result_stat4$comparisons
"IL9Tx_IL9R_AQ - IL9Tx_IL9R_PR" "IL9Tx_IL9R_AQ - IL9Tx_IL9R_WT" "IL9Tx_IL9R_PR - IL9Tx_IL9R_WT"


#pstat5

stat5a_values <- GetAssayData(merged_obj_subset, assay = "regulonAUC", slot = "data")["Stat5a(+)",]
kruskal_result_stat5a <- kruskal.test(stat5a_values ~ Idents(merged_obj_subset))
dunn_result_stat5a <- dunn.test(stat5a_values, g = Idents(merged_obj_subset), method = "holm")
#get the stat values
dunn_result_stat5a$P
2.133819e-178 1.505131e-128  5.324533e-47

dunn_result_stat5a$comparisons
"IL9Tx_IL9R_AQ - IL9Tx_IL9R_PR" "IL9Tx_IL9R_AQ - IL9Tx_IL9R_WT" "IL9Tx_IL9R_PR - IL9Tx_IL9R_WT"







#ggsave("~/Desktop/Figure1-IL9/stats_pyscenic.pdf", width = 28, height = 14, units = "cm")

#List of STAT regulons
stat_regulons <- c("Stat1(+)", "Stat2(+)", "Stat3(+)", "Stat4(+)", "Stat5a(+)", "Stat5b(+)","Stat6(+)")

#Get the data in the correct form
#Extract the AUC data
#it is stored in regulon assay in counts and data slot
#head(merged_obj@assays$regulonAUC$counts@x)

auc_data <- as.data.frame(t(merged_obj_subset@assays$regulonAUC$data[stat_regulons,]))

#Add treatment information
auc_data$Treatment <- Idents(merged_obj_subset)

#Reshape the data for plotting
auc_data_long <- auc_data %>%
  pivot_longer(cols = all_of(stat_regulons), names_to = "Regulon", values_to = "AUC")

#PLOTTING

#Now perform the analysis on subset object 
cluster_colors_treatment <- c("IL9Tx_IL9R_WT" = "blue",
                              "IL9Tx_IL9R_PR" = "#F39B7FFF",
                              "IL9Tx_IL9R_AQ" = "#91D1C2FF")
#Create the ridge plot
ggplot(auc_data_long, aes(x = AUC, y = Regulon, fill = Treatment)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9, rel_min_height = 0.01) +
  facet_wrap(~ Treatment, ncol = 1) +
  scale_fill_manual(values = cluster_colors_treatment) +
  theme_ridges() +
  labs(title = "STAT Regulon Activity Across Treatments",
       x = "AUC", y = "STAT Regulon") +
  theme(legend.position = "none",  # Remove legend as it's redundant with facets
        strip.text = element_text(face = "bold"),  # Make facet labels bold
        axis.title = element_text(face = "bold"),  # Make axis titles bold
        plot.title = element_text(face = "bold", hjust = 0.5))  # Center and bold the main title


#Save top regulons from scenic by cluster names
# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(readr)

#Load the seurat object
merged_obj <- readRDS("~/Desktop/Kayla_scData/merged_obj.rds")
merged_obj <- JoinLayers(merged_obj)
DefaultAssay(merged_obj) <-"RNA"

#by treatment
Idents(merged_obj) <-"seurat_clusters"

#merged_obj <- RenameIdents(merged_obj, `0` = "Tem/rm", `1` = "T early activated", `2` = "Teff-1/prolif",
#                           `3` = "Teff-2/prolif", `4` = "Teff-3/pex/naïve(mix)", `5` = "Teff-4/naïve/em(mix)", `6` = "Teff-5(Stat1/Gzmb/Prf1)", `7` = "Tscm-like(Tcf7)",`8`="Tmyeloid-like",`9`="Tprolif")
merged_obj$seurat_clusters <- Idents(merged_obj) 
#new regulons script
regulons <- read_csv("~/Desktop/Kayla_scData/scenic/regulons_new.csv")
head(regulons)

#Load the AUC Matrix
regulonAUC <- importAUCfromText("~/Desktop/Kayla_scData/scenic/auc_mtx_new.csv")  ##New 
Idents(merged_obj) <- "seurat_clusters"
cellInfo <- data.frame(seurat_clusters=Idents(merged_obj))
cellInfo <- as.data.frame(cellInfo)

cellsPerCluster <- split(rownames(cellInfo), cellInfo$seurat_clusters) 

#cellsPerCluster <- as.data.frame(cellsPerCluster)
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = F, scale=T))

temp <- as.data.frame(regulonActivity_byCellType)
temp$tf <- rownames(temp) #This is the relative exrichment score This is achieved using AUCell, which provides an AUC value indicating the relative expression of the genes in the regulon (as a whole) in each cell.

#####
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "seurat_clusters", "RelativeActivity")
topRegulators$seurat_clusters <- factor(as.character(topRegulators$seurat_clusters))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators) #These are the top regulators by treatment group 
write.csv(topRegulators, "~/Desktop/Kayla_scData/scenic/allRegulators_by_clusters_noPruning.csv", row.names = FALSE)

topRegulators_by_clusters <- topRegulators %>%
  group_by(seurat_clusters) %>%
  top_n(n = 250, wt = RelativeActivity)

write.csv(topRegulators_by_clusters, "~/Desktop/Kayla_scData/scenic/topRegulators_by_clusters_noPruning.csv", row.names = FALSE)

