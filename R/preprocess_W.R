#' @title Align network matrix rows/columns with gene list from DE results
#'
#' @description Ensures that the network matrix includes all genes present in the DE results.
#' Missing genes are added as zero rows/columns. The final matrix is subset and
#' reordered to match the gene order in `de`.
#'
#' @param de Data frame or matrix with gene names as row names.
#' @param W Square matrix with row and column names (gene identifiers).
#'
#' @return A square matrix with rows/columns named after all genes in `de`,
#'         ordered identically to `rownames(de)`.
#'
#' @examples
#' \dontrun{
#' de <- data.frame(pvalue = runif(10), row.names = paste0("gene", 1:10))
#' # Create network matrix for 8 of these genes (2 missing)
#' W <- matrix(runif(64), 8, 8, 
#'             dimnames = list(paste0("gene", 1:8), paste0("gene", 1:8)))
#' # Align: adds gene9 and gene10 as zero rows/columns
#' W_aligned <- preprocess_W(de, W)
#' }
#' @export
preprocess_W <- function(de, W){
  W_zeros <- rownames(de)[!rownames(de) %in% rownames(W)]
  W1 <- Matrix(0, nrow = dim(W)[1]+length(W_zeros), ncol = dim(W)[1]+length(W_zeros), sparse = TRUE)
  W1[1:dim(W)[1],1:dim(W)[1]] <- W
  rownames(W1) <- c(rownames(W), W_zeros)
  colnames(W1) <- c(rownames(W), W_zeros)
  W1 <- as.matrix(W1)

  
  W1 <- W1[rownames(de),rownames(de)]
}
  
