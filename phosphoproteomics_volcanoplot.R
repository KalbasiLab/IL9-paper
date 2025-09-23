
#Phosphoproteomic data 
library(readr)
library(EnhancedVolcano)
#read in the data 
#IL9R_IL9_vs_IL9R_Mock
IL9R_IL9_vs_IL9R_Mock <- read_csv("Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL9_vs_IL9R_Mock.csv")
#UNIPROT, GENE_SYMBOL, phosposition, phosAAposition, average_intensity, logFC, P, adj_P

IL9R_IL9_vs_IL9R_Mock_df<- as.data.frame(IL9R_IL9_vs_IL9R_Mock[c("uniprot","description","gene_symbol","Gene.names","P.Value","adj.P.Val","logFC","gene_symbol")])
#write.csv(IL9R_IL9_vs_IL9R_Mock_df, "~/Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL9_vs_IL9R_Mock_re.csv", quote = F)

#remove the special character from the gene_symbol col
IL9R_IL9_vs_IL9R_Mock_df$gene_symbol <-gsub("[[:punct:]]", "", IL9R_IL9_vs_IL9R_Mock_df$gene_symbol)
#create unique rownames because many gene names are present at positions in the gene symbol 
rownames(IL9R_IL9_vs_IL9R_Mock_df) = make.names(IL9R_IL9_vs_IL9R_Mock_df$gene_symbol, unique=TRUE) 


#plotting
keyvals <- rep('black', nrow(IL9R_IL9_vs_IL9R_Mock_df))
names(keyvals) <- rep('Not significant', nrow(IL9R_IL9_vs_IL9R_Mock_df))

# Identify significant points
significant <- which(IL9R_IL9_vs_IL9R_Mock_df$adj.P.Val < 0.05 & abs(IL9R_IL9_vs_IL9R_Mock_df$logFC) > 0.5)
keyvals[significant] <- 'darkred'
names(keyvals)[significant] <- 'Significant'

# Make column names unique
names(IL9R_IL9_vs_IL9R_Mock_df) <- make.names(names(IL9R_IL9_vs_IL9R_Mock_df), unique = TRUE)


#Create the plot
EnhancedVolcano(IL9R_IL9_vs_IL9R_Mock_df,
                lab = IL9R_IL9_vs_IL9R_Mock_df$gene_symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'IL9R_IL9 vs IL9R_Mock',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                #labSize = 6.0,  # Reduced label size to fit more labels
                colCustom = keyvals,
                colAlpha = 1.5,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                ylim = c(0,3), 
                xlim = c(-3,7))
                #selectLab = IL9R_IL9_vs_IL9R_Mock_df$gene_symbol[significant],
                #max.overlaps = Inf,  # Allow unlimited label overlaps
                #boxedLabels = TRUE,  # Add boxes around labels to improve readability
                #min.segment.length = 0,  # Allow shorter connector lines
                #maxoverlapsConnectors = Inf,
                #cutoffLineType = 'blank')  # Allow unlimited connector overlaps

ggsave("~/Desktop/Figure1-IL9/Figure3/phosphoproteomics_IL9vsMock_withoutLabel.pdf", width = 18, height = 18, units = "cm")

#revision plot #add the gene positions labels 
# For manual labeling of a specific gene
# Start with the normal gene labels for all points
labels <- IL9R_IL9_vs_IL9R_Mock_df$gene_symbol

# Find the Stat4 entries
specific_gene <- "Stat4"
specific_idx <- which(IL9R_IL9_vs_IL9R_Mock_df$gene_symbol == specific_gene)

# Get their log2FC values
fc_values <- IL9R_IL9_vs_IL9R_Mock_df$logFC[specific_idx]

# Customize only the Stat4 labels with their log2FC values
idx_fc_6_25 <- specific_idx[which(abs(fc_values - 6.25) == min(abs(fc_values - 6.25)))]
idx_fc_1_2 <- specific_idx[which(abs(fc_values - 1.2) == min(abs(fc_values - 1.2)))]

labels[idx_fc_6_25] <- paste0(specific_gene, " (Y694)")
labels[idx_fc_1_2] <- paste0(specific_gene, " (S722)")

# Create the plot with all gene labels, including the custom Stat4 labels
EnhancedVolcano(IL9R_IL9_vs_IL9R_Mock_df,
                lab = labels,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'IL9R_IL9 vs IL9R_Mock',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                #labSize = 3.5,
                #boxedLabels = TRUE,
                drawConnectors = TRUE,
                #max.overlaps = 10,
                colCustom = keyvals,
                colAlpha = 1.5,
                legendPosition = 'right',
                #legendLabSize = 10,
                widthConnectors = 0.25,
                ylim = c(0,3), 
                xlim = c(-3,7))

ggsave("~/Desktop/Figure1-IL9/Figure3/phosphoproteomics_IL9vsMock_withoutLabel_updated.pdf", width = 17, height = 16.5, units = "cm")


########################second###################
##IL9R_IL2_vs_IL9R_Mock
IL9R_IL2_vs_IL9R_Mock <- read_csv("Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL2_vs_IL9R_Mock.csv")
head(IL9R_IL2_vs_IL9R_Mock)

IL9R_IL2_vs_IL9R_Mock_df<- as.data.frame(IL9R_IL2_vs_IL9R_Mock[c("uniprot","description","gene_symbol","Gene.names","P.Value","adj.P.Val","logFC","gene_symbol")])

write.csv(IL9R_IL2_vs_IL9R_Mock_df, "~/Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL2_vs_IL9R_Mock_re.csv", quote = F)


#IL9R_IL2_vs_IL9R_Mock_df<- as.data.frame(IL9R_IL2_vs_IL9R_Mock[c("adj.P.Val","logFC","gene_symbol")])
#remove the special character from the gene_symbol col

IL9R_IL2_vs_IL9R_Mock_df$gene_symbol <-gsub("[[:punct:]]", "", IL9R_IL2_vs_IL9R_Mock_df$gene_symbol)
#create unique rownames because many gene names are present at positions in the gene symbol 
rownames(IL9R_IL2_vs_IL9R_Mock_df) = make.names(IL9R_IL2_vs_IL9R_Mock_df$gene_symbol, unique=TRUE) 

#plotting
#####Only get the significant ones
keyvals <- rep('black', nrow(IL9R_IL2_vs_IL9R_Mock_df))
names(keyvals) <- rep('Not significant', nrow(IL9R_IL2_vs_IL9R_Mock_df))

# Identify significant points
significant <- which(IL9R_IL2_vs_IL9R_Mock_df$adj.P.Val < 0.05 & abs(IL9R_IL2_vs_IL9R_Mock_df$logFC) > 0.5)
#228
keyvals[significant] <- 'darkred'

names(keyvals)[significant] <- 'Significant'

#Create the plot
selectLab <- c("Stat5a")
EnhancedVolcano(IL9R_IL2_vs_IL9R_Mock_df,
                lab = IL9R_IL2_vs_IL9R_Mock_df$gene_symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'IL9R_IL2 vs IL9R_Mock',
                pCutoff = 0.6,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 6.0,  # Reduced label size to fit more labels
                colCustom = keyvals,
                colAlpha = 1.5,
                legendPosition = 'right',
                selectLab = selectLab,
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                ylim = c(0,7), 
                xlim = c(-3,7))

ggsave("~/Desktop/Figure1-IL9/Figure3/phosphoproteomics_IL2vsMock.pdf", width = 28, height = 30, units = "cm")

########################Third###################
#IL9R_IL9_vs_IL9R_IL2
IL9R_IL9_vs_IL9R_IL2 <- read_csv("Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL9_vs_IL9R_IL2.csv")
head(IL9R_IL9_vs_IL9R_IL2)

IL9R_IL9_vs_IL9R_IL2_df<- as.data.frame(IL9R_IL9_vs_IL9R_IL2[c("uniprot","description","gene_symbol","Gene.names","P.Value","adj.P.Val","logFC","gene_symbol")])

write.csv(IL9R_IL9_vs_IL9R_IL2_df, "~/Desktop/Phosphoproteomic/KalbasiA_032124/phospho_analysis/02_DA_results/per_contrast_results/IL9R_IL9_vs_IL9R_IL2_re.csv", quote = F)

#IL9R_IL9_vs_IL9R_IL2_df<- as.data.frame(IL9R_IL9_vs_IL9R_IL2[c("adj.P.Val","logFC","gene_symbol")])
#remove the special character from the gene_symbol col

IL9R_IL9_vs_IL9R_IL2_df$gene_symbol <-gsub("[[:punct:]]", "", IL9R_IL9_vs_IL9R_IL2_df$gene_symbol)
#create unique rownames because many gene names are present at positions in the gene symbol 
rownames(IL9R_IL9_vs_IL9R_IL2_df) = make.names(IL9R_IL9_vs_IL9R_IL2_df$gene_symbol, unique=TRUE) 

#Plotting
keyvals <- rep('black', nrow(IL9R_IL9_vs_IL9R_IL2_df))
names(keyvals) <- rep('Not significant', nrow(IL9R_IL9_vs_IL9R_IL2_df))

# Identify significant points
significant <- which(IL9R_IL9_vs_IL9R_IL2_df$adj.P.Val < 0.05 & abs(IL9R_IL9_vs_IL9R_IL2_df$logFC) > 0.5)
#264 
keyvals[significant] <- 'darkred'
names(keyvals)[significant] <- 'Significant'

#Create the plot
EnhancedVolcano(IL9R_IL9_vs_IL9R_IL2_df,
                lab = IL9R_IL9_vs_IL9R_IL2_df$gene_symbol,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'IL9R_IL9 vs IL9R_IL2',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 6.0,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                ylim = c(0,5),
                xlim = c(-3, 6))
                #selectLab = rownames(IL9R_IL9_vs_IL9R_IL2_df)[significant],
                #cutoffLineType = 'blank',
                #max.overlaps = Inf,  # Allow unlimited label overlaps
                #boxedLabels = TRUE,  # Add boxes around labels to improve readability
                #min.segment.length = 0,  # Allow shorter connector lines
                #maxoverlapsConnectors = Inf)

ggsave("~/Desktop/Figure1-IL9/Figure3/phosphoproteomics_IL9vsIL2_withoutLabel.pdf", width = 18, height = 18, units = "cm")



