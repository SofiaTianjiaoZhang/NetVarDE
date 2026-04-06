#' @title NetVarDE_Infer with Gene-Specific Network Sensitivity (βᵢ)
#'
#' @description Implements network-enhanced differential expression analysis with
#' gene-specific network sensitivity parameters (βᵢ). This allows different
#' genes to have different sensitivities to network influence.
#'
#' @param p_values Numeric vector of p-values (length n).
#' @param network_matrix Symmetric adjacency matrix (n x n) with non-negative weights.
#' @param gene_names Optional character vector of gene names (default NULL).
#' @param method Method for estimating βᵢ:
#'   \itemize{
#'     \item \code{"variational"}: full variational inference (most accurate but slow).
#'     \item \code{"covariate_variational"}: uses external covariate to inform βᵢ prior.
#'   }
#' @param beta_prior_mean Prior mean for βᵢ (default: 1.0).
#' @param beta_prior_var Prior variance for βᵢ (default: 0.5).
#' @param covariate Optional covariate vector for \code{"covariate_variational"} method.
#' @param max_iter Maximum EM iterations (default: 100).
#' @param tol Convergence tolerance (default: 1e-6).
#' @param estimate_mixture_method Method for estimating null/alternative distributions:
#'   \code{"beta"} or \code{"kde"} (default: \code{"beta"}).
#' @param verbose Print progress messages (default: TRUE).
#' @param seed Random seed for reproducibility (default: 123).
#'
#' @return A list of class \code{nelfdr_beta_i_result} with components:
#'   \item{gene}{Gene names (if provided).}
#'   \item{p_value}{Original p-values.}
#'   \item{posterior}{Posterior probability of being differentially expressed.}
#'   \item{local_pvalue}{1 - posterior.}
#'   \item{parameters}{List of estimated parameters (alpha, beta_i, etc.).}
#'   \item{beta_i}{Gene-specific network sensitivity parameters.}
#'   \item{beta_i_stats}{Summary statistics of βᵢ (mean, sd, min, max, CV).}
#'   \item{convergence}{Convergence information (iterations, converged flag, final delta).}
#'   \item{settings}{List of model settings used.}
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n_genes <- 500
#' p_vals <- c(rbeta(100, 0.1, 1), runif(400))
#' W <- matrix(rbinom(n_genes^2, 1, 0.02), n_genes, n_genes)
#' W <- (W + t(W)) / 2
#' diag(W) <- 0
#'
#' result <- NetVarDE_Infer(p_vals, W, method = "variational")
#'
#' }
#' @export
NetVarDE_Infer <- function(p_values,
                                      network_matrix,
                                      gene_names=NULL,
                                      method = "variational",
                                      beta_prior_mean = 1.0,
                                      beta_prior_var = 0.5,
                                      covariate=c(),
                                      max_iter = 100,
                                      tol = 1e-6,
                                      estimate_mixture_method="beta",
                                      verbose = TRUE,
                                      seed = 123) {
  
  set.seed(seed)
  
  # ========== 1. Input Validation ==========
  if (verbose) cat("Step 1: Validating inputs...\n")
  
  if (!is.numeric(p_values)) stop("p_values must be numeric")
  if (any(p_values < 0 | p_values > 1)) {
    stop("p_values must be in [0, 1]")
  }
  
  n_genes <- length(p_values)
  
  if (!is.matrix(network_matrix)) {
    stop("network_matrix must be a matrix")
  }
  if (nrow(network_matrix) != ncol(network_matrix)) {
    stop("network_matrix must be square")
  }
  if (nrow(network_matrix) != n_genes) {
    stop("Dimensions of p_values and network_matrix don't match")
  }
  
  valid_methods <- c("variational","covariate_variational")
  if (!method %in% valid_methods) {
    stop(paste("method must be one of:", paste(valid_methods, collapse = ", ")))
  }
  if (method == "covariate_variational") {
    if (length(covariate) != n_genes) stop("covariate length must equal number of genes")
  }
  
  # Normalize network matrix
  W <- normalize_network_matrix(network_matrix, verbose = verbose, method = "row")
  
  # Handle zero p-values
  p_adj <- p_values
  p_adj[p_adj == 0] <- min(.Machine$double.eps, min(p_adj[p_adj > 0]*0.1, na.rm = TRUE))
  p_adj[p_adj == 1] <- max(1 - .Machine$double.eps, max(p_adj[p_adj < 1], na.rm = TRUE))
  
  
  # ========== 2. Estimate Mixture Components ==========
  if (verbose) cat("Step 2: Estimating null and alternative distributions...\n")
  
  mixture <- estimate_mixture_components(p_adj, method = estimate_mixture_method, verbose = FALSE)
  
  # ========== 3. Estimate Gene-Specific βᵢ ==========
  if (verbose) cat(sprintf("Step 3: Estimating gene-specific βᵢ using '%s' method...\n", method))
  
  beta_i <- switch(method,
                   variational = estimate_beta_variational(p_adj, W, mixture,
                                                           beta_prior_mean,beta_prior_var,
                                                           max_iter, tol, verbose),
                   covariate_variational = estimate_beta_covariate_variational(p_adj, W, mixture,covariate, 
                                                                             beta_prior_var, max_iter, tol)
  )
  
  # ========== 4. Run Inference with Gene-Specific βᵢ ==========
  if (verbose) cat("Step 4: Running inference with gene-specific βᵢ...\n")
  
  result <- ising_mean_field_beta_i_stable(p_adj, W, mixture, beta_i,
                                           max_iter, tol, verbose)
  
  # ========== 5. Compute Additional Diagnostics ==========
  if (verbose) cat("Step 5: Computing diagnostics...\n")
  
  # Compute local pvalue
  local_pvalue <- 1 - result$posterior
  
  # ========== 6. Prepare Output ==========
  output <- list(
    # Core results
    gene = gene_names,
    p_value = p_values,
    posterior = result$posterior,
    local_pvalue = local_pvalue,
    
    # Gene-specific parameters
    parameters = result$parameters,
    
    # Beta_i information
    beta_i = beta_i,
    beta_i_stats = list(
      mean = mean(beta_i),
      sd = sd(beta_i),
      min = min(beta_i),
      max = max(beta_i),
      cv = sd(beta_i) / mean(beta_i)  # Coefficient of variation
    ),
    
    # Convergence info
    convergence = result$convergence,
    
    # Model settings
    settings = list(
      method = method,
      beta_prior_mean = beta_prior_mean,
      beta_prior_var = beta_prior_var,
      n_genes = n_genes,
      seed = seed
    )
  )
  
  class(output) <- "nelfdr_beta_i_result"
  
  if (verbose) {
    cat("\n========================================\n")
    cat("Gene-Specific βᵢ Analysis Complete\n")
    cat("========================================\n")
    cat(sprintf("Method: %s\n", method))
    cat(sprintf("βᵢ statistics:\n"))
    cat(sprintf("  Mean: %.3f, SD: %.3f\n", mean(beta_i), sd(beta_i)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(beta_i), max(beta_i)))
    cat(sprintf("  CV: %.3f\n", sd(beta_i)/mean(beta_i)))
    if (!is.null(result$convergence)) {
      cat(sprintf("Converged in %d iterations\n", result$convergence$iterations))
    }
    cat("========================================\n")
  }
  
  return(output)
}