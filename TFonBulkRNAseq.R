
#Identifying tfs from bulk rna seq data 
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

#decoupleR analysis
net <- get_collectri(organism='mouse', split_complexes=FALSE)
net

#Compare the aq vs Wt
res_WTvsAQ <- results(dds1, contrast=c("condition","IL9RPosIL.9","IL9R.AQPosIL.9")) #--Comparison1

res_WTvsAQ$gene <- rownames(res_WTvsAQ)
res_WTvsAQ <- as.data.frame(res_WTvsAQ)

#Do DEG between the AQvs wt and PR vs WT
rownames(res_WTvsAQ) <- make.unique(toupper(rownames(res_WTvsAQ)))
res_WTvsAQ$gene <-toupper(res_WTvsAQ$gene)
#res_il9vsil2$gene <-toupper(res_il9vsil2$gene)
deg <- res_WTvsAQ %>%
  as.data.frame() %>%
  dplyr::select(gene, log2FoldChange, stat, pvalue) %>% 
  dplyr::filter(!is.na(stat)) 

#Set rownames and convert to matrix
deg_matrix <- as.matrix(deg[, -1]) # Exclude the first column for matrix conversion
rownames(deg_matrix) <- deg$gene

# Check the result
head(deg_matrix)

contrast_acts <- run_ulm(mat=deg_matrix[, 'stat', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

############
n_tfs <- 30
# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  dplyr::filter(source %in% tfs)

#plot
ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score, fill = score > 0)) + 
  geom_bar(stat="identity", width=0.8, color="black", size=0.1) +
  scale_fill_manual(values=c("darkblue", "darkred"), labels = c("AQ","WT")) +
  coord_flip() +
  labs(x = "Transcription factors", 
       y = "Enrichment Score", 
       title = "TFs between WT vs AQ",
       fill = "Group") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,color="black"),
        axis.title = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        axis.text.x = element_text(size = 12,color="black"),
        legend.title = element_text(size = 12,color="black"),
        legend.text = element_text(size = 12,color="black"))

ggsave("~/Desktop/Figure1-IL9/Figure3/TF_WTvsAQ_top30.pdf", width = 16, height = 18, units = "cm")



