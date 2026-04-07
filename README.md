## Quick Start

Here is an example of using `NetVarDE` to analyze differential expression with a gene network.

```r
# Load package
library(NetVarDE)

# ------------------------------------------------------------
# 1. Prepare DE results 
# ------------------------------------------------------------
# For demonstration, we provided an example DE data in data folder
# In practice, you would load your own DE results.
DE_example <- read.csv("data/DE_example.csv", header = TRUE, row.names = 1)

# (Optional) Remove mitochondria-related genes if working with snRNA-seq data
MT <- read.csv("data/Human.MitoCarta3.0.csv")
DE_example <- DE_example[!rownames(DE_example) %in% MT$Symbol,]

# ------------------------------------------------------------
# 2. Prepare a gene network (e.g., PPI network from STRING database)
# ------------------------------------------------------------
l <- load("data/Prior_Matrix.Rdata")

# ------------------------------------------------------------
# 3. Align network matrix with DE genes
# ------------------------------------------------------------
W_aligned <- preprocess_W(DE_example, W)
dim(W_aligned)  # should match number of genes in 'DE_example'

# ------------------------------------------------------------
# 4. (Optional) Sparse network selection via power-law fit
# ------------------------------------------------------------
sparse_result <- power_law_selection(W_aligned,
                                     min_density = 0.01,
                                     max_density = 0.3,
                                     n_bins = 20,
                                     Rsquared_cutoff = 0.85)
W_sparse <- sparse_result$W_sparse

# ------------------------------------------------------------
# 5. Run network‑enhanced DE inference
# ------------------------------------------------------------
results <- NetVarDE_Infer(p_values = DE_example$p_val,
                          network_matrix = W_sparse,
                          gene = rownames(DE_example))

# Create summary data frame
results_df <- data.frame(gene = results$gene,
                         p_value = results$p_value,
                         NetVarDE_local_pvalue = results$local_pvalue)

# View top genes with strongest evidence
head(results_df[order(results_df$NetVarDE_local_pvalue), ])
