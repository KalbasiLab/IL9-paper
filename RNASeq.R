#Bulk RNA-seq 
#August 5th 
library(tidyverse)
library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("ggplot2")
library(EnhancedVolcano)
library(ggfortify)
library(PCAtools)
library(ggfortify)

###lOAD THE counts data and create the metadata to do the DEG
counts_data <- read.delim('~/Desktop/Kayla-RNASeq/gene_counts.tsv')
#counts_data <- read.delim('~/Desktop/Kayla-RNASeq/salmon.merged.gene_counts_length_scaled.tsv')
head(counts_data)
#Load the counts data as tximport

#rownames(counts_data) <- counts_data$gene_name
rownames(counts_data) <- make.names(counts_data$gene_name, unique = TRUE)
countdata <- counts_data[,-(1:2)]
#rownames(countdata) <- counts_data[,2]
head(countdata)
colnames(countdata)

#create the metadata
treatment <- c("IL9R.AQPosIL.9_1","IL9R.AQPosIL.9_2","IL9R.AQPosIL.9_3","IL9R.PRPosIL.9_1", "IL9R.PRPosIL.9_2","IL9R.PRPosIL.9_3",
               "IL9R.py5XPosIL.9_1","IL9R.py5XPosIL.9_2","IL9R.py5XPosIL.9_3","IL9RPosIL.2_1","IL9RPosIL.2_2","IL9RPosIL.2_3",
               "IL9RPosIL.9_1", "IL9RPosIL.9_2","IL9RPosIL.9_3","o9RPosMSA.oIL2_1","o9RPosMSA.oIL2_2","o9RPosMSA.oIL2_3") 

replicates <- rep(c(1, 2,3), times = 6)

condition <- c("IL9R.AQPosIL.9","IL9R.AQPosIL.9" , "IL9R.AQPosIL.9" , "IL9R.PRPosIL.9" , "IL9R.PRPosIL.9" , "IL9R.PRPosIL.9" ,"IL9R.py5XPosIL.9","IL9R.py5XPosIL.9","IL9R.py5XPosIL.9",
               "IL9RPosIL.2" ,  "IL9RPosIL.2"  , "IL9RPosIL.2"   ,"IL9RPosIL.9"  ,"IL9RPosIL.9" , "IL9RPosIL.9", "o9RPosMSA.oIL2" ,"o9RPosMSA.oIL2" , "o9RPosMSA.oIL2")

metadata <- data.frame(
  treatment = treatment,
  condition = condition
)

table(colnames(countdata)==metadata$treatment)

#Step 2: construct a DESeqDataSet object ----------
dds <- DESeqDataSetFromMatrix(countData = round(countdata),
                              colData = metadata,
                              design = ~ condition)

#Filter out the low expressed genes  below 10 reads per samples 
#this is considered the minimal pre-filtering applied before DESeq2 analysis to keep only rows that have at least 10 reads total)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds1 <- dds[keep,]

#DESeq2
dds1 <- DESeq(dds1)
head(assay(dds1)) 
res_table <- results(dds1)
resultsNames(dds1)
#summary(res_table)

#Compare the Il9 vs Il2  
res_il9vsil2 <- results(dds1, contrast=c("condition","IL9RPosIL.9","IL9RPosIL.2")) #--Comparison1
res_il9vsil2$gene <- rownames(res_il9vsil2)
res_il9vsil2 <- as.data.frame(res_il9vsil2)

resOrdered <- res_il9vsil2[order(res_il9vsil2$pvalue),]
summary(resOrdered)


#Change the color and formatting
# Create a custom color vector
keyvals <- rep('black', nrow(res_il9vsil2))
names(keyvals) <- rep('Not significant', nrow(res_il9vsil2))

# Identify significant points
significant <- which(res_il9vsil2$pvalue < 10e-6 & abs(res_il9vsil2$log2FoldChange) > 2)
keyvals[significant] <- 'darkred'
names(keyvals)[significant] <- 'Significant'

#Create the plot
EnhancedVolcano(res_il9vsil2,
                #lab = IL9R_IL9_vs_IL9R_IL2_df$gene_symbol,
                lab=NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'IL9R_IL9 vs IL9R_IL2',
                pCutoff = 10e-6,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 6.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                #selectLab = rownames(IL9R_IL9_vs_IL9R_IL2_df)[significant],
                cutoffLineType = 'blank')


ggsave("~/Desktop/Figure1-IL9/Figure3/il9vsil2.pdf", width = 20, height = 16, units = "cm")


##Compare the orthoIL9 vs Il9
res_il9rvso9r<- results(dds1, contrast=c("condition","IL9RPosIL.9","o9RPosMSA.oIL2")) #--Comparision2
res_il9rvso9r$gene <- rownames(res_il9rvso9r)
res_il9rvso9r <- as.data.frame(res_il9rvso9r)
resOrdered <- res_il9rvso9r[order(res_il9rvso9r$pvalue),]
summary(resOrdered)
#Plotting
EnhancedVolcano(res_il9rvso9r,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'IL9vsoIL9',
                pointSize = 2.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'darkred'),
                colAlpha = 1)

ggsave("~/Desktop/Figure1-IL9/Figure3/il9vsoil9.pdf", width = 14, height = 18, units = "cm")


# Create a custom color vector
keyvals <- rep('black', nrow(res_il9rvso9r))
names(keyvals) <- rep('Not significant', nrow(res_il9rvso9r))

#Identify significant points
significant <- which(res_il9rvso9r$pvalue < 10e-6 & abs(res_il9rvso9r$log2FoldChange) > 2)
keyvals[significant] <- 'darkred'
names(keyvals)[significant] <- 'Significant'

#Create the plot
EnhancedVolcano(res_il9rvso9r,
                #lab = IL9R_IL9_vs_IL9R_IL2_df$gene_symbol,
                lab=NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'IL9 vs oIL9',
                pCutoff = 10e-6,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 6.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                #selectLab = rownames(IL9R_IL9_vs_IL9R_IL2_df)[significant],
                cutoffLineType = 'blank')

ggsave("~/Desktop/Figure1-IL9/Figure3/il9vsoil9.pdf", width = 20, height = 16, units = "cm")




#Ortho9r vs IL2
res_o9Rvsil2 <- results(dds1, contrast=c("condition","o9RPosMSA.oIL2","IL9RPosIL.2")) #--Comparision3
res_o9Rvsil2$gene <- rownames(res_o9Rvsil2)
res_o9Rvsil2 <- as.data.frame(res_o9Rvsil2)
#write.csv(res_o9Rvsil2, "~/Desktop/Kayla-RNASeq/DEG/res_o9Rvsil2.csv", quote = F)

resOrdered <- res_o9Rvsil2[order(res_o9Rvsil2$pvalue),]
summary(resOrdered)

#Change the color and formatting
# Create a custom color vector
keyvals <- rep('black', nrow(res_il9vsil2))
names(keyvals) <- rep('Not significant', nrow(res_il9vsil2))

#Identify significant points
significant <- which(res_o9Rvsil2$pvalue < 10e-6 & abs(res_o9Rvsil2$log2FoldChange) > 2)
keyvals[significant] <- 'darkred'
names(keyvals)[significant] <- 'Significant'

#Create the plot
EnhancedVolcano(res_o9Rvsil2,
                #lab = res_o9Rvsil2$gene_symbol,
                lab=NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'o9R vs IL9R_IL2',
                pCutoff = 10e-6,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 6.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                #selectLab = rownames(IL9R_IL9_vs_IL9R_IL2_df)[significant],
                cutoffLineType = 'blank')

ggsave("~/Desktop/Figure1-IL9/Figure3/oil9vsil2.pdf", width = 20, height = 16, units = "cm")













#Get all the significant genes and then do the overlapping #this gene list would include the both the upregulated and downregulated genes -- but only the significant genes 
res_o9Rvsil2_sig <- res_o9Rvsil2 %>%
  dplyr::filter(padj < 0.05)
#7166
#res_o9rvsIl9_sig <- res_o9rvsIl9 %>%
#  dplyr::filter(padj < 0.05)
#3031
res_il9rvso9R_sig <- res_il9rvso9r %>%
  dplyr::filter(padj < 0.05)

res_il9vsil2_sig <- res_il9vsil2 %>%
  dplyr::filter(padj < 0.05)
#7673

#The venn diagram is going to be between the genes that are significant and upregulated
##Get the venn diagram for the DEG between conditions
library(dplyr)
library(ggVennDiagram)
library("ggvenn") 

# Extract gene lists
o9Rvsil2 <- res_o9Rvsil2_sig$gene
il9vso9R<- res_o9rvsIl9_sig$gene
il9vsil2 <- res_il9vsil2_sig$gene

# Create a list of gene sets
gene_sets <- list(
  o9Rvsil2 = o9Rvsil2,
  il9vso9R = il9vso9R,
  il9vsil2 = il9vsil2)

ggvenn(gene_sets,show_percentage=FALSE)


####
##Venn diagram of genes upregulated by IL9RvsIL2 and genes upregulated by o9RvIL2
#Il9vsil2 (IL9RPosIL VS IL9RPosIL.2 --res_il9vsil2)
#o9Rvsil2  (o9RPosMSA.oIL2 VS IL9RPosIL.2 --res_o9Rvsil2)

#Get only the upregulated genes  and significant 
log2fc_threshold <- 1.5
padj_threshold <- 0.05
# Subset the data to include only positive and significant genes
res_o9Rvsil2_pos_sig <- subset(res_o9Rvsil2_sig,padj < padj_threshold & log2FoldChange > log2fc_threshold)
res_o9Rvsil2_pos_sig

# Subset the data to include only positive and significant genes
res_il9vsil2_pos_sig <- subset(res_il9vsil2_sig,padj < padj_threshold & log2FoldChange > log2fc_threshold)
res_il9vsil2_pos_sig

#Extract gene lists
o9Rvsil2 <- res_o9Rvsil2_pos_sig$gene
il9vsil2 <- res_il9vsil2_pos_sig$gene

# Create a list of gene sets
gene_sets <- list(
  o9Rvsil2 = o9Rvsil2,
  il9vsil2 = il9vsil2)

ggvenn(gene_sets,show_percentage=FALSE)



###Genes common between IL9R (IL9RPosIL) vs o9R (o9RPosMSA.oIL2) but different in IL2 (IL9RPosIL.2) group? 
#res_o9Rvsil2$gene #DEG between o9RPosMSA.oIL2 vs IL9RPosIL.2
#res_o9Rvsil9$gene #DEG between o9RPosMSA.oIL2 vs IL9RPosIL
#res_il9vsil2$gene #Deg between IL9RPosIL vs IL9RPosIL.2

common_genes_o9Rvsil2_and_il9vsil2 <- intersect(res_o9Rvsil9$gene, res_il9vsil2$gene)
#2390
# Find the genes that are different from res_il9vsil2$gene
result_genes_il9ANDorthoVSIL2 <- setdiff(common_genes_o9Rvsil2_and_il9vsil2, res_il9vsil2$gene)

length(result_genes_il9ANDorthoVSIL2)
#641 genes are different between DEG (orthoIL9 VS IL9)  and (IL9R VS IL2) group but

#Question2 
#thresholds for significance and positive log2FC
padj_threshold <- 0.05
log2fc_threshold <- 0

# Subset the data to include only positive and significant genes
res_il9vsil2_pos_sig <- subset(res_il9vsil2,padj < padj_threshold & log2FoldChange > log2fc_threshold)
res_il9vsil2_pos_sig
#3698

#Now get the genes that are downregulated in 
#padj_threshold <- 0.05
log2fc_threshold <- 0

#o9r vs IL2
#res_o9Rvsil2
# Subset the data to include only negative (less than 0 log2FC ) and significant genes 
#res_o9Rvsil2_NEGATIVE_sig <- subset(res_o9Rvsil2, padj < padj_threshold & log2FoldChange < log2fc_threshold)
res_o9Rvsil2_NEGATIVE_sig <- subset(res_o9Rvsil2, padj < padj_threshold &log2FoldChange < log2fc_threshold)
dim(res_o9Rvsil2_NEGATIVE_sig)
#3635

#Intersect to get the genes that are going up in IL9VSIL2 group but down in o9RvsIL2 group
common_genesIL9VSIL2_o9RvsIL2 <- intersect(res_il9vsil2_pos_sig$gene, res_o9Rvsil2_NEGATIVE_sig$gene)
common_genesIL9VSIL2_o9RvsIL2
#66

#Heatmap of these genes across three conditions
dds1_vst = vst(dds1, blind=FALSE) #vst normalised counts data from DESEQ2 object 
normalised_counts_data<- assay(dds1_vst)
normalised_counts_data<- as.data.frame(normalised_counts_data) #this is my normalised counts 
normalised_counts_data_subset <- normalised_counts_data[, !colnames(normalised_counts_data) %in% c("IL9R.py5XPosIL.9_1", "IL9R.py5XPosIL.9_2","IL9R.py5XPosIL.9_3","IL9R.AQPosIL.9_1","IL9R.AQPosIL.9_2","IL9R.AQPosIL.9_3",
                                                                                                   "IL9R.PRPosIL.9_1","IL9R.PRPosIL.9_2","IL9R.PRPosIL.9_3")]
#BIOCARTA_IL12 heatmap

# Specify the desired column order
desired_order <- c(
  "IL9RPosIL.2_1", "IL9RPosIL.2_2", "IL9RPosIL.2_3", 
  "o9RPosMSA.oIL2_1", "o9RPosMSA.oIL2_2", "o9RPosMSA.oIL2_3",
  "IL9RPosIL.9_1", "IL9RPosIL.9_2", "IL9RPosIL.9_3"
)

# Reorder the columns in TActivation_fate matrix
normalised_counts_data <- normalised_counts_data[, desired_order]

genes_BIOCARTA_IL12<- c("Cd247","Tyk2","Jak2","Il8r1","Mapk14","Cd3e","Ifng","Il12rb1","Il12rb2","Il18","Cd3g","Cxcr3","Ccr5","Jun","Stat4","Map2k6","Cd3d","Mapk8","Tyk2","Etv5") #another gene is doubled

BIOCARTA_IL12 <- normalised_counts_data[rownames(normalised_counts_data) %in% genes_BIOCARTA_IL12, ]
library(pheatmap)

P1 <- pheatmap(as.matrix(BIOCARTA_IL12), 
               scale = "row",  # Scale rows (genes)
               clustering_distance_rows = "euclidean",  # Clustering distance metric for rows
               clustering_distance_cols = "euclidean",
               cluster_cols=TRUE,# Clustering distance metric for columns
               cluster_rows=TRUE,
               clustering_method = "complete",  # Clustering method
               main = "DEG Heatmap BIOCARTA_IL12 ",  # Title of the heatmap
               fontsize = 12,  # Font size
               cellwidth = 20,  # Cell width
               cellheight = 20,color=colorRampPalette(c("navy", "white", "red"))(50))
P1



