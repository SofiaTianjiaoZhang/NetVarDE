#' @title Ising Mean Field Inference with Gene-Specific βᵢ
#'
#' @description Performs mean-field variational inference for the Ising model
#' with gene-specific network sensitivity parameters βᵢ. This function estimates
#' the posterior probability of differential expression for each gene given
#' the p-values, network structure, and estimated βᵢ.
#'
#' @param p_values Numeric vector of p-values (length n).
#' @param W Normalized adjacency matrix (n x n). Should be row-normalized
#'   (each row sums to 1) or symmetric normalized. Use \code{normalize_network_matrix}
#'   to preprocess.
#' @param mixture List containing mixture components (f0, f1, pi0, pi1) from
#'   \code{estimate_mixture_components}.
#' @param beta_i Numeric vector of gene-specific βᵢ (length n), typically
#'   estimated by \code{estimate_beta_variational} or similar.
#' @param max_iter Maximum number of EM iterations (default 100).
#' @param tol Convergence tolerance based on change in posterior probabilities
#'   (default 1e-6).
#' @param verbose If TRUE, print progress messages every 10 iterations.
#'
#' @return A list with components:
#'   \item{posterior}{Posterior probability of differential expression (vector).}
#'   \item{parameters}{List containing alpha, beta_i, beta_global, beta_sd, pi0.}
#'   \item{convergence}{List with iterations, converged flag, final delta.}
#'
#' @keywords internal
#' 
#' 
ising_mean_field_beta_i_stable <- function(p_values, W, mixture, beta_i,
                                           max_iter = 100, tol = 1e-6,
                                           verbose = TRUE) {
  
  n <- length(p_values)
  
  if (verbose) cat("  Initializing inference with gene-specific βᵢ...\n")
  
  # Initialize posterior probabilities from local FDR
  lfdr_init <- compute_local_fdr_initial(p_values, mixture)
  gamma <- 1 - lfdr_init
  
  # Clamp initial values to avoid extremes
  gamma <- pmin(pmax(gamma, 0.001), 0.999)
  
  alpha <- 0
  
  # Use the provided normalized network matrix
  
  W_norm <- W
  
  # Precompute log likelihood ratios and clamp to stable range
  log_lr <- compute_log_likelihood_ratio(p_values, mixture)
  log_lr <- pmin(pmax(log_lr, -10), 10)  # Clamp to [-10, 10]
  
  
  for (iter in 1:max_iter) {
    gamma_old <- gamma
    
    # --- E-step: Update posterior probabilities ---
    neighbor_influence <- W_norm %*% gamma
    
    # Numerically stable calculation of eta
    eta <- log_lr + alpha + beta_i * neighbor_influence
    
    # Numerically stable sigmoid with piecewise evaluation
    gamma_new <- numeric(n)
    
    # For large positive eta
    large_pos <- eta > 100
    gamma_new[large_pos] <- 1 - exp(-eta[large_pos])
    
    # For large negative eta
    large_neg <- eta < -100
    gamma_new[large_neg] <- exp(eta[large_neg])
    
    # For moderate eta
    moderate <- !(large_pos | large_neg)
    gamma_new[moderate] <- 1 / (1 + exp(-eta[moderate]))
    
    # Ensure numerical stability
    gamma_new <- pmin(pmax(gamma_new, 1e-10), 1 - 1e-10)
    gamma <- gamma_new
    
    # --- M-step: Update global parameters ---
    # Update alpha using stable Newton-Raphson (max 5 iterations)
    alpha <- update_alpha_stable(gamma, p_values, mixture, beta_i, W,
                                 max_iter = 5, tol = 1e-4)
    
    # Clamp alpha to avoid extreme values
    alpha <- pmin(pmax(alpha, -10), 10)
    
    # --- Convergence check ---
    delta <- max(abs(gamma - gamma_old))
    
    if (verbose && iter %% 10 == 0) {
      cat(sprintf("    Iter %3d: delta = %.2e, alpha = %.3f\n",
                  iter, delta, alpha))
    }
    
    if (delta < tol) {
      if (verbose) cat(sprintf("  Converged after %d iterations\n", iter))
      break
    }
  }
  
  ## Prepare convergence info
  convergence <- list(
    iterations = iter,
    converged = iter < max_iter,
    final_delta = delta
  )
  
  parameters <- list(
    alpha = alpha,
    beta_i = beta_i,
    beta_global = mean(beta_i),
    beta_sd = sd(beta_i),
    pi0 = mixture$pi0
  )
  
  return(list(
    posterior = gamma,
    parameters = parameters,
    convergence = convergence
  ))
}