
#Correlation plot

#Load  the data 
counts_data <- read.delim('~/Desktop/Kayla-RNASeq/gene_counts.tsv')
head(counts_data)
#make the samples as colnames
rownames(counts_data) <- make.names(counts_data$gene_name, unique = TRUE)
countdata <- counts_data[,-(1:2)]
#rownames(countdata) <- counts_data[,2]
head(countdata)
colnames(countdata)

#For the metadata
#create the metadata
treatment <- c("IL9R.AQPosIL.9_1","IL9R.AQPosIL.9_2","IL9R.AQPosIL.9_3","IL9R.PRPosIL.9_1",
               "IL9R.PRPosIL.9_2","IL9R.PRPosIL.9_3",
               "IL9R.py5XPosIL.9_1","IL9R.py5XPosIL.9_2","IL9R.py5XPosIL.9_3","IL9RPosIL.2_1","IL9RPosIL.2_2","IL9RPosIL.2_3",
               "IL9RPosIL.9_1", "IL9RPosIL.9_2","IL9RPosIL.9_3","o9RPosMSA.oIL2_1","o9RPosMSA.oIL2_2","o9RPosMSA.oIL2_3") 

replicates <- rep(c(1, 2,3), times = 6)

condition <- c("IL9R.AQPosIL.9","IL9R.AQPosIL.9" , "IL9R.AQPosIL.9" , "IL9R.PRPosIL.9" , "IL9R.PRPosIL.9" , "IL9R.PRPosIL.9" ,"IL9R.py5XPosIL.9","IL9R.py5XPosIL.9","IL9R.py5XPosIL.9",
               "IL9RPosIL.2" ,  "IL9RPosIL.2"  , "IL9RPosIL.2"   ,"IL9RPosIL"  ,"IL9RPosIL" , "IL9RPosIL", "o9RPosMSA.oIL2" ,"o9RPosMSA.oIL2" , "o9RPosMSA.oIL2")

metadata <- data.frame(
  treatment = treatment,
  condition = condition
)

#write.csv(metadata, "~/Desktop/Kayla-RNASeq/metadata.csv", row.names = FALSE)
table(colnames(countdata)==metadata$treatment)

# Step 1: Add pseudocount (usually 1)
#counts_plus1 <- counts_data + 1

# Step 2: Log2 transform
#log_normalized_counts <- log2(counts_plus1)
#final_data <- as.data.frame(log_normalized_counts)
#Load the excel file 
library(readxl) # For reading Excel files
# Read the Excel file
processed_pSTATdata <- read_excel("~/Desktop/Kayla-RNASeq/phosphoseq-link/2024-10-03-pSTAT-Correlation.xlsx")
# Convert to a dataframe
processed_pSTATdata <- as.data.frame(processed_pSTATdata)
rownames(processed_pSTATdata) <- processed_pSTATdata$Emax
#processed_pSTATdata$...1 <- NULL
processed_pSTATdata <- processed_pSTATdata[,-(1)]
#Rename the columns according to the RNA-seq counts data to stack them

name_mapping <- c(
  "WT + IL9" = "IL9RPosIL.9_1",
  "...5" = "IL9RPosIL.9_2",
  "AQ" = "IL9R.AQPosIL.9_1",
  "...7" = "IL9R.AQPosIL.9_2",
  "PR" = "IL9R.PRPosIL.9_1",
  "...9" = "IL9R.PRPosIL.9_2",
  "5x" = "IL9R.py5XPosIL.9_1",
  "...13" = "IL9R.py5XPosIL.9_2",
  "o9R" = "o9RPosMSA.oIL2_1",
  "...11" = "o9RPosMSA.oIL2_2",
  "WT + IL2" = "IL9RPosIL.2_1",
  "...3" = "IL9RPosIL.2_2"
)

#rename the colnames to match the RNA-seq 
colnames(processed_pSTATdata) <- name_mapping[colnames(processed_pSTATdata)]

# Start with the original data
result <- processed_pSTATdata

# Define the groups
#groups <- c("IL9R.AQPosIL.9", "IL9R.PRPosIL.9", "IL9RPosIL.2", "IL9RPosIL.9", "IL9R.py5XPosIL.9", "o9RPosMSA.oIL2")
# Create a list of column pairs
column_pairs <- list(
  c("IL9RPosIL.2_1", "IL9RPosIL.2_2"),
  c("IL9RPosIL.9_1", "IL9RPosIL.9_2"),
  c("IL9R.AQPosIL.9_1", "IL9R.AQPosIL.9_2"),
  c("IL9R.PRPosIL.9_1", "IL9R.PRPosIL.9_2"),
  c("o9RPosMSA.oIL2_1", "o9RPosMSA.oIL2_2"),
  c("IL9R.py5XPosIL.9_1", "IL9R.py5XPosIL.9_2")
)

# Function to calculate mean of two columns
calculate_mean <- function(df, col1, col2) {
  rowMeans(df[, c(col1, col2)], na.rm = TRUE)
}

#Calculate means for each pair and add to the dataframe
for (pair in column_pairs) {
  new_col_name <- gsub("_1", "_mean", pair[1])
  processed_pSTATdata[[new_col_name]] <- calculate_mean(processed_pSTATdata, pair[1], pair[2])
}

# Remove the original paired columns
columns_to_remove <- unlist(column_pairs)
processed_pSTATdata <- processed_pSTATdata[, !(names(processed_pSTATdata) %in% columns_to_remove)]

# View the resulting dataframe
head(processed_pSTATdata)

#triplicate 
# Get the current mean columns
mean_columns <- grep("_mean$", names(processed_pSTATdata), value = TRUE)

# Create a new dataframe to store the triplicates
new_data <- data.frame(row.names = rownames(processed_pSTATdata))

#For each mean column, create three identical columns
for (col in mean_columns) {
  base_name <- sub("_mean$", "", col)
  for (i in 1:3) {
    new_col_name <- paste0(base_name, "_", i)
    new_data[[new_col_name]] <- processed_pSTATdata[[col]]
  }
}

#Replace the original dataframe with the new one
processed_pSTATdata <- new_data

#View the resulting dataframe
head(processed_pSTATdata)

#Define the desired column order
desired_order <- c(
  "IL9R.AQPosIL.9_1", "IL9R.AQPosIL.9_2", "IL9R.AQPosIL.9_3",
  "IL9R.PRPosIL.9_1", "IL9R.PRPosIL.9_2", "IL9R.PRPosIL.9_3",
  "IL9R.py5XPosIL.9_1", "IL9R.py5XPosIL.9_2", "IL9R.py5XPosIL.9_3",
  "IL9RPosIL.2_1", "IL9RPosIL.2_2", "IL9RPosIL.2_3",
  "IL9RPosIL.9_1", "IL9RPosIL.9_2", "IL9RPosIL.9_3",
  "o9RPosMSA.oIL2_1", "o9RPosMSA.oIL2_2", "o9RPosMSA.oIL2_3"
)

# Reorder the columns
processed_pSTATdata <- processed_pSTATdata[, desired_order]
# View the result
head(processed_pSTATdata)

#Combine these two datasets
common_columns <- intersect(colnames(countdata), colnames(processed_pSTATdata))
colnames(countdata) == colnames(processed_pSTATdata)
#1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# Log normalize phospho-STAT data
#processed_pSTATdata_plus1 <- processed_pSTATdata + 1
#log_normalized_pSTATdata <- log2(processed_pSTATdata_plus1)
#log normalise these two together 


#Now just take the average of the of each pstat values so that is consistent  
# Define the groups
#this is the mean value of pstat across the samples (like all AQs would  have the same value, (AQ1+AQ2+AQ3)/3 AND likewise for all other samples) 
#Combine by rows (concatenate)
combined_Data<- rbind(countdata, processed_pSTATdata)
expression_matrix <- as.data.frame(combined_Data)
genes_of_interest <- rownames(expression_matrix)

#write.csv(expression_matrix, "~/Desktop/Kayla-RNASeq/phosphoseq-link/expression_matrix.csv", quote = F)

#Log normalize phospho-STAT and RNA-seq together log2 normalisation is most common in gene expression data 
combined_Data_plus1 <- combined_Data + 1  #log of zero is undefined so handling 0 values
log_normalized_combineddata <- log2(combined_Data_plus1)

####Correlation######
pSTAT1_expression <- expression_matrix["pSTAT1", ] #this is the rownmaes and this could be the gene that we would be using to correlate other genes. 
genes_expression <- expression_matrix[genes_of_interest, ] 

#Extract expression data for the specific gene (here it could be pstats) and the genes of interest
selected_genes <- c("pSTAT1", genes_of_interest) #Get the PSTAT1 as the first col in the expression matrix 
subset_expression_matrix <- expression_matrix[selected_genes, ] #This would contain the expression matrix of all the genes including an additional row at the top whihc corresponts to any of the pstats

#Transpose the expression matrix to switch rows and columns
transposed_expression_matrix <- t(subset_expression_matrix) #Divide this into three dataframe and transpose them independently  

#Now correlate the pstat1 which is present as the columns with other columns values (genes)
# Initialize an empty vector to store the correlation values
correlation_values <- c()
# Loop through all genes (excluding Gene1 column which is pstat1) to calculate the correlation with Gene1
for (gene in colnames(transposed_expression_matrix)[-1]) {
  correlation <- cor(transposed_expression_matrix[, "pSTAT1"], transposed_expression_matrix[, gene], use = "complete.obs")
  correlation_values <- c(correlation_values, correlation)
}

#Create a dataframe to store the correlations
correlation_df <- data.frame(
  Gene = colnames(transposed_expression_matrix)[-1], #The first col contains the pstat1 (which is what we want to correlate every gene with)
  Correlation_with_pSTAT1= correlation_values
)

#Print the resulting dataframe
head(correlation_df) #This contains the correlation of pstat1 with all the genes in the dataframe
#write.csv(correlation_df, "~/Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT1.csv", quote = F)
#Load the correlation csv for plotting
library(readr)
correlation_df_pSTAT1 <- read_csv("Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT1.csv")
#expression_matrix <- read_csv("Desktop/Kayla-RNASeq/phosphoseq-link/expression_matrix.csv")

#Order the correlations in descending order and get the top 100 genes with are correlated with the pstats (here pstat1)
top_100_genes <- correlation_df_pSTAT1[order(-correlation_df_pSTAT1$Correlation_with_pSTAT1), ][1:100, ]

#Print the top 100 genes with their correlation values with pSTAT1
top_100_genes <- head(top_100_genes, 100)
print(top_100_genes)

pSTAT1_expr <- expression_matrix["pSTAT1", ]
top_100_expr <- expression_matrix[top_100_genes$Gene, ]

#Replace "pSTAT1.1" with "pSTAT1" in the top_100_genes dataframe
top_100_genes$Gene <- ifelse(top_100_genes$Gene == "pSTAT1.1", "pSTAT1", top_100_genes$Gene)

#Verify the change
print(top_100_genes[top_100_genes$Gene == "pSTAT1", ])
top_100_expr <- expression_matrix[top_100_genes$Gene, ] #Expression matrix should have genes as rownames

#get the module score for these 100 genes which is basically calculated by taking the means/ average of the columns(samples) 
#so that every sample has one value for expression of all 100 genes. note: always drop the pstat gene expression. 
module_expr <- colMeans(top_100_expr[top_100_genes$Gene != "pSTAT1", ])

print(head(module_expr))

plot_data <- data.frame(
  Sample = colnames(expression_matrix),
  pSTAT1_expression = as.numeric(expression_matrix["pSTAT1", ]),
  Module_expression = module_expr
)

#Plot this
ggplot(plot_data, aes(x = pSTAT1_expression, y = Module_expression)) +
  geom_point(aes(color = Sample), size = 3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(
    x = "pSTAT1 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT1 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )+
  scale_color_manual(values = c(
    "IL9R.AQPosIL.9_1" = "#91D1C2FF", "IL9R.AQPosIL.9_2" = "#91D1C2FF", "IL9R.AQPosIL.9_3" = "#91D1C2FF",
    "IL9R.PRPosIL.9_1" = "#F39B7FFF", "IL9R.PRPosIL.9_2" = "#F39B7FFF", "IL9R.PRPosIL.9_3" = "#F39B7FFF",
    "IL9R.py5XPosIL.9_1" = "#B9DBF4", "IL9R.py5XPosIL.9_2" = "#B9DBF4", "IL9R.py5XPosIL.9_3" = "#B9DBF4",
    "IL9RPosIL.2_1" = "black", "IL9RPosIL.2_2" = "black", "IL9RPosIL.2_3" = "black",
    "IL9RPosIL.9_1" = "blue", "IL9RPosIL.9_2" = "blue", "IL9RPosIL.9_3" = "blue",
    "o9RPosMSA.oIL2_1" = "#4DBBD5FF", "o9RPosMSA.oIL2_2" = "#4DBBD5FF", "o9RPosMSA.oIL2_3" = "#4DBBD5FF"
  ))

#Calculate R-squared
lm_model <- lm(Module_expression ~ pSTAT1_expression, data = plot_data)
r_squared <- summary(lm_model)$r.squared

# Create the plot with R-squared annotation

##Final plotting --Nove 4th 
# First, let's create a simplified version of the Sample column
plot_data$Sample_Group <- gsub("_[1-3]$", "", plot_data$Sample)

#Now, let's create the plot
ggplot(plot_data, aes(x = pSTAT1_expression, y = Module_expression)) +
  geom_point(aes(color = Sample_Group), size = 37) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dotted", size = 2) +
  theme_minimal() +
  labs(
    x = "pSTAT1 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT1 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(axis.title.x = element_text(color = "black", size = 60),
        axis.text.x = element_text(hjust = 1, color = "black", size = 90),
        axis.text.y = element_text(lineheight = 0, color = "black", size = 90),
        axis.title.y = element_text(color = "black", size = 60),
        text = element_text(size = 60),
        legend.title = element_text(color = "black", size = 60),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -1, size = 60, vjust=2.5),
        axis.title = element_text(size = 60),
        legend.position = "right"
  ) +
  annotate(
    "text",
    x = max(plot_data$pSTAT1_expression),
    y = min(plot_data$Module_expression),
    label = sprintf("R² = %.2f", r_squared),
    hjust = 1,
    vjust = 0,
    size = 50
  ) +
  scale_color_manual(
    values = c(
      "IL9R.AQPosIL.9" = "#91D1C2FF",
      "IL9R.PRPosIL.9" = "#F39B7FFF",
      "IL9R.py5XPosIL.9" = "#B9DBF4",
      "IL9RPosIL.2" = "black",
      "IL9RPosIL.9" = "blue",
      "o9RPosMSA.oIL2" = "#4DBBD5FF"
    ),
    name = "Sample Group"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))

ggsave("~/Desktop/Figure1-IL9/pstat1_correlation_new.pdf", width = 78, height = 80, units = "cm")






#Pstat3
correlation_df_pSTAT3 <- read_csv("Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT3.csv")

#Order the correlations in descending order and get the top 100 genes with are correlated with the pstats (here pstat1)
top_100_genes_pstat3 <- correlation_df_pSTAT3[order(-correlation_df_pSTAT3$Correlation_with_pSTAT3), ][1:100, ]

pSTAT3_expr <- expression_matrix["pSTAT3", ]
top_100_expr <- expression_matrix[top_100_genes_pstat3$Gene, ]

# Replace "pSTAT3.1" with "pSTAT3" in the top_100_genes dataframe
top_100_genes_pstat3$Gene <- ifelse(top_100_genes_pstat3$Gene == "pSTAT3.1", "pSTAT3", top_100_genes_pstat3$Gene)

# Verify the change
print(top_100_genes_pstat3[top_100_genes_pstat3$Gene == "pSTAT3", ])

top_100_expr <- expression_matrix[top_100_genes_pstat3$Gene, ]

#get the module score for these 100 genes which is basically calculated by taking the means/ average of the columns(samples) 
#so that every sample has one value for expression of all 100 genes 

module_expr <- colMeans(top_100_expr[top_100_genes_pstat3$Gene != "pSTAT3", ])

# Check the result
print(head(module_expr))


plot_data <- data.frame(
  Sample = colnames(expression_matrix),
  pSTAT3_expression = as.numeric(expression_matrix["pSTAT3", ]),
  Module_expression = module_expr
)


#Plot this
ggplot(plot_data, aes(x = pSTAT3_expression, y = Module_expression)) +
  geom_point(aes(color = Sample), size = 3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(
    x = "pSTAT3 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT1 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )+
  scale_color_manual(values = c(
    "IL9R.AQPosIL.9_1" = "#91D1C2FF", "IL9R.AQPosIL.9_2" = "#91D1C2FF", "IL9R.AQPosIL.9_3" = "#91D1C2FF",
    "IL9R.PRPosIL.9_1" = "#F39B7FFF", "IL9R.PRPosIL.9_2" = "#F39B7FFF", "IL9R.PRPosIL.9_3" = "#F39B7FFF",
    "IL9R.py5XPosIL.9_1" = "#B9DBF4", "IL9R.py5XPosIL.9_2" = "#B9DBF4", "IL9R.py5XPosIL.9_3" = "#B9DBF4",
    "IL9RPosIL.2_1" = "black", "IL9RPosIL.2_2" = "black", "IL9RPosIL.2_3" = "black",
    "IL9RPosIL.9_1" = "blue", "IL9RPosIL.9_2" = "blue", "IL9RPosIL.9_3" = "blue",
    "o9RPosMSA.oIL2_1" = "#4DBBD5FF", "o9RPosMSA.oIL2_2" = "#4DBBD5FF", "o9RPosMSA.oIL2_3" = "#4DBBD5FF"
  ))

#Calculate R-squared
lm_model <- lm(Module_expression ~ pSTAT3_expression, data = plot_data)
r_squared <- summary(lm_model)$r.squared

#Now, let's create the plot
ggplot(plot_data, aes(x = pSTAT3_expression, y = Module_expression)) +
  geom_point(aes(color = Sample_Group), size = 37) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dotted", size = 2) +
  theme_minimal() +
  labs(
    x = "pSTAT3 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT3 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(axis.title.x = element_text(color = "black", size = 60),
        axis.text.x = element_text(hjust = 1, color = "black", size = 90),
        axis.text.y = element_text(lineheight = 0, color = "black", size = 90),
        axis.title.y = element_text(color = "black", size = 60),
        text = element_text(size = 60),
        legend.title = element_text(color = "black", size = 60),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -1, size = 60, vjust=2.5),
        axis.title = element_text(size = 60),
        legend.position = "right"
  ) +
  annotate(
    "text",
    x = max(plot_data$pSTAT3_expression),
    y = min(plot_data$Module_expression),
    label = sprintf("R² = %.2f", r_squared),
    hjust = 1,
    vjust = 0,
    size = 50
  ) +
  scale_color_manual(
    values = c(
      "IL9R.AQPosIL.9" = "#91D1C2FF",
      "IL9R.PRPosIL.9" = "#F39B7FFF",
      "IL9R.py5XPosIL.9" = "#B9DBF4",
      "IL9RPosIL.2" = "black",
      "IL9RPosIL.9" = "blue",
      "o9RPosMSA.oIL2" = "#4DBBD5FF"
    ),
    name = "Sample Group"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))


ggsave("~/Desktop/Figure1-IL9/pstat3_correlation_new.pdf", width = 78, height = 80, units = "cm")


#pSTAT4
correlation_df_pSTAT4 <- read_csv("Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT4.csv")
#Order the correlations in descending order and get the top 100 genes with are correlated with the pstats (here pstat1)
top_100_genes <- correlation_df_pSTAT4[order(-correlation_df_pSTAT4$Correlation_with_pSTAT4), ][1:100, ]

#Print the top 100 genes with their correlation values with pSTAT1
top_100_genes <- head(top_100_genes, 100)
head(top_100_genes)

#correlation
pSTAT4_expr <- expression_matrix["pSTAT4", ]
top_100_expr <- expression_matrix[top_100_genes$Gene, ]

# Replace "pSTAT4.1" with "pSTAT4" in the top_100_genes dataframe
top_100_genes$Gene <- ifelse(top_100_genes$Gene == "pSTAT4.1", "pSTAT4", top_100_genes$Gene)

# Verify the change
print(top_100_genes[top_100_genes$Gene == "pSTAT4", ])

top_100_expr <- expression_matrix[top_100_genes$Gene, ]

#get the module score for these 100 genes which is basically calculated by taking the means/ average of the columns(samples) 
#so that every sample has one value for expression of all 100 genes 
module_expr <- colMeans(top_100_expr[top_100_genes$Gene != "pSTAT4", ])

# Check the result
print(head(module_expr))

plot_data <- data.frame(
  Sample = colnames(expression_matrix),
  pSTAT4_expression = as.numeric(expression_matrix["pSTAT4", ]),
  Module_expression = module_expr
)

#Plot this
ggplot(plot_data, aes(x = pSTAT4_expression, y = Module_expression)) +
  geom_point(aes(color = Sample), size = 3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(
    x = "pSTAT1 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT4 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )+
  scale_color_manual(values = c(
    "IL9R.AQPosIL.9_1" = "#91D1C2FF", "IL9R.AQPosIL.9_2" = "#91D1C2FF", "IL9R.AQPosIL.9_3" = "#91D1C2FF",
    "IL9R.PRPosIL.9_1" = "#F39B7FFF", "IL9R.PRPosIL.9_2" = "#F39B7FFF", "IL9R.PRPosIL.9_3" = "#F39B7FFF",
    "IL9R.py5XPosIL.9_1" = "#B9DBF4", "IL9R.py5XPosIL.9_2" = "#B9DBF4", "IL9R.py5XPosIL.9_3" = "#B9DBF4",
    "IL9RPosIL.2_1" = "black", "IL9RPosIL.2_2" = "black", "IL9RPosIL.2_3" = "black",
    "IL9RPosIL.9_1" = "blue", "IL9RPosIL.9_2" = "blue", "IL9RPosIL.9_3" = "blue",
    "o9RPosMSA.oIL2_1" = "#4DBBD5FF", "o9RPosMSA.oIL2_2" = "#4DBBD5FF", "o9RPosMSA.oIL2_3" = "#4DBBD5FF"
  ))

#Calculate R-squared
lm_model <- lm(Module_expression ~ pSTAT4_expression, data = plot_data)
r_squared <- summary(lm_model)$r.squared


#Plot stat4
plot_data$Sample_Group <- gsub("_[1-3]$", "", plot_data$Sample)

#Now, let's create the plot
ggplot(plot_data, aes(x = pSTAT4_expression, y = Module_expression)) +
  geom_point(aes(color = Sample_Group), size = 37) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dotted", size = 2) +
  theme_minimal() +
  labs(
    x = "pSTAT4 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT4 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(axis.title.x = element_text(color = "black", size = 60),
        axis.text.x = element_text(hjust = 1, color = "black", size = 90),
        axis.text.y = element_text(lineheight = 0, color = "black", size = 90),
        axis.title.y = element_text(color = "black", size = 60),
        text = element_text(size = 60),
        legend.title = element_text(color = "black", size = 60),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -1, size = 60, vjust=2.5),
        axis.title = element_text(size = 60),
        legend.position = "right"
  ) +
  annotate(
    "text",
    x = max(plot_data$pSTAT4_expression),
    y = min(plot_data$Module_expression),
    label = sprintf("R² = %.2f", r_squared),
    hjust = 1,
    vjust = 0,
    size = 50
  ) +
  scale_color_manual(
    values = c(
      "IL9R.AQPosIL.9" = "#91D1C2FF",
      "IL9R.PRPosIL.9" = "#F39B7FFF",
      "IL9R.py5XPosIL.9" = "#B9DBF4",
      "IL9RPosIL.2" = "black",
      "IL9RPosIL.9" = "blue",
      "o9RPosMSA.oIL2" = "#4DBBD5FF"
    ),
    name = "Sample Group"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))


ggsave("~/Desktop/Figure1-IL9/pstat4_correlation_new.pdf", width = 78, height = 80, units = "cm")



#pSTAT5
#Print the top 100 genes with their correlation values with pSTAT4
pSTAT5_expression <- expression_matrix["pSTAT5", ] #this is the rownmaes and this could be the gene that we would be using to correlate other genes. 
genes_expression <- expression_matrix[genes_of_interest, ] 

#Extract expression data for the specific gene (here it could be pstats) and the genes of interest
selected_genes <- c("pSTAT5", genes_of_interest) #Get the PSTAT5 as the first col in the expression matrix 
subset_expression_matrix <- expression_matrix[selected_genes, ] #This would contain the expression matrix of all the genes including an additional row at the top whihc corresponts to any of the pstats

#Transpose the expression matrix to switch rows and columns
transposed_expression_matrix <- t(subset_expression_matrix) 
#Now correlate the pstat1 which is present as the columns with other columns values (genes)
# Initialize an empty vector to store the correlation values
correlation_values_pstat5 <- c()
# Loop through all genes (excluding Gene1 which is pstat5) to calculate the correlation with Gene1
for (gene in colnames(transposed_expression_matrix)[-1]) {
  correlation <- cor(transposed_expression_matrix[, "pSTAT5"], transposed_expression_matrix[, gene], use = "complete.obs")
  correlation_values_pstat5 <- c(correlation_values_pstat5, correlation)
}

# Create a dataframe to store the correlations
correlation_df_pstat5 <- data.frame(
  Gene = colnames(transposed_expression_matrix)[-1], #The first col contains the pstat1 (which is what we want to correlate every gene with)
  Correlation_with_pSTAT5= correlation_values_pstat5
)

#write.csv(correlation_df_pstat5, "~/Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT5.csv", quote = F)

correlation_df_pSTAT5 <- read_csv("Desktop/Kayla-RNASeq/phosphoseq-link/correlation_df_pSTAT5.csv")

#Order the correlations in descending order and get the top 100 genes with are correlated with the pstats (here pstat1)
top_100_genes <- correlation_df_pSTAT5[order(-correlation_df_pSTAT5$Correlation_with_pSTAT5), ][1:100, ]
#Module score of these top 100 correlated genes
pSTAT5_expr <- expression_matrix["pSTAT5", ]
top_100_expr <- expression_matrix[top_100_genes$Gene, ]

# Replace "pSTAT5.1" with "pSTAT5" in the top_100_genes dataframe
top_100_genes$Gene <- ifelse(top_100_genes$Gene == "pSTAT5.1", "pSTAT5", top_100_genes$Gene)

# Verify the change
print(top_100_genes[top_100_genes$Gene == "pSTAT5", ])

top_100_expr <- expression_matrix[top_100_genes$Gene, ]

#get the module score for these 100 genes which is basically calculated by taking the means/ average of the columns(samples) 
#so that every sample has one value for expression of all 100 genes. note: always drop the pstat gene expression. 
module_expr <- colMeans(top_100_expr[top_100_genes$Gene != "pSTAT5", ])

print(head(module_expr))

plot_data <- data.frame(
  Sample = colnames(expression_matrix),
  pSTAT5_expression = as.numeric(expression_matrix["pSTAT5", ]),
  Module_expression = module_expr
)


#Plot this
ggplot(plot_data, aes(x = pSTAT5_expression, y = Module_expression)) +
  geom_point(aes(color = Sample), size = 3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(
    x = "pSTAT5 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT1 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(
    axis.title.x = element_text(color = "black", size = 60),
    axis.text.x = element_text(hjust = 1, color = "black", size =60),
    axis.text.y = element_text(lineheight = 0, color = "black", size = 60),
    axis.title.y = element_text(color = "black", size = 60),
    text = element_text(size = 60),
    legend.title = element_text(color = "black", size = 60),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size=60),
    axis.title = element_text(size=60),
    legend.position = "right"
  )+
  scale_color_manual(values = c(
    "IL9R.AQPosIL.9_1" = "#91D1C2FF", "IL9R.AQPosIL.9_2" = "#91D1C2FF", "IL9R.AQPosIL.9_3" = "#91D1C2FF",
    "IL9R.PRPosIL.9_1" = "#F39B7FFF", "IL9R.PRPosIL.9_2" = "#F39B7FFF", "IL9R.PRPosIL.9_3" = "#F39B7FFF",
    "IL9R.py5XPosIL.9_1" = "#B9DBF4", "IL9R.py5XPosIL.9_2" = "#B9DBF4", "IL9R.py5XPosIL.9_3" = "#B9DBF4",
    "IL9RPosIL.2_1" = "black", "IL9RPosIL.2_2" = "black", "IL9RPosIL.2_3" = "black",
    "IL9RPosIL.9_1" = "blue", "IL9RPosIL.9_2" = "blue", "IL9RPosIL.9_3" = "blue",
    "o9RPosMSA.oIL2_1" = "#4DBBD5FF", "o9RPosMSA.oIL2_2" = "#4DBBD5FF", "o9RPosMSA.oIL2_3" = "#4DBBD5FF"
  ))

#Calculate R-squared
lm_model <- lm(Module_expression ~ pSTAT5_expression, data = plot_data)
r_squared <- summary(lm_model)$r.squared

#Plot stat4
plot_data$Sample_Group <- gsub("_[1-3]$", "", plot_data$Sample)

#Now, let's create the plot
ggplot(plot_data, aes(x = pSTAT5_expression, y = Module_expression)) +
  geom_point(aes(color = Sample_Group), size = 37) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dotted", size = 2) +
  theme_minimal() +
  labs(
    x = "pSTAT5 Expression",
    y = "Top 100 Correlated Genes Module Expression",
    title = "pSTAT5 Expression vs. Top 100 Correlated Genes Module"
  ) +
  theme(axis.title.x = element_text(color = "black", size = 60),
        axis.text.x = element_text(hjust = 1, color = "black", size = 90),
        axis.text.y = element_text(lineheight = 0, color = "black", size = 90),
        axis.title.y = element_text(color = "black", size = 60),
        text = element_text(size = 60),
        legend.title = element_text(color = "black", size = 60),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -1, size = 60, vjust=2.5),
        axis.title = element_text(size = 60),
        legend.position = "right"
  ) +
  annotate(
    "text",
    x = max(plot_data$pSTAT5_expression),
    y = min(plot_data$Module_expression),
    label = sprintf("R² = %.2f", r_squared),
    hjust = 1,
    vjust = 0,
    size = 50
  ) +
  scale_color_manual(
    values = c(
      "IL9R.AQPosIL.9" = "#91D1C2FF",
      "IL9R.PRPosIL.9" = "#F39B7FFF",
      "IL9R.py5XPosIL.9" = "#B9DBF4",
      "IL9RPosIL.2" = "black",
      "IL9RPosIL.9" = "blue",
      "o9RPosMSA.oIL2" = "#4DBBD5FF"
    ),
    name = "Sample Group"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1))


ggsave("~/Desktop/Figure1-IL9/pstat5_correlation_new.pdf", width = 78, height = 80, units = "cm")
