#' @title Covariate-Informed Variational EM for βᵢ
#'
#' @description Estimates gene-specific network sensitivity parameters βᵢ using
#' variational inference with a hierarchical prior that depends on an external
#' covariate. The prior mean for each βᵢ is modeled as a linear function of the
#' covariate, allowing borrowing of information across genes.
#'
#' @param p_values Numeric vector of p-values (length n).
#' @param W Normalized adjacency matrix (n x n).
#' @param mixture List containing mixture components (f0, f1, pi0, pi1) from
#'   \code{estimate_mixture_components}.
#' @param covariate_x Numeric vector of external covariate (e.g., gene sparsity,
#'   expression level) with length n. Used to inform the prior mean of βᵢ.
#' @param prior_var Variance for the βᵢ hierarchical layer (default 0.5).
#' @param max_iter Maximum number of EM iterations (default 50).
#' @param tol Convergence tolerance based on change in gamma (default 1e-4).
#'
#' @return Numeric vector of estimated βᵢ (length n).
#'
#' @keywords internal
#' 
#' 
estimate_beta_covariate_variational <- function(p_values, W, mixture, covariate_x, 
                                                prior_var = 0.5, 
                                                max_iter = 50, tol = 1e-4) {
  n <- length(p_values)
  
  # --- 1. Pre-computation and Initialization ---
  # Normalize network weights once to save computation
  W_norm <- W 
  W_norm_sq <- (W)^2
  
  # Initial gene states (gamma) and log-likelihood ratios
  gamma <- 1 - compute_local_fdr_initial(p_values, mixture)
  log_lr <- compute_log_likelihood_ratio(p_values, mixture)
  alpha <- 0 # Global baseline intercept
  
  # Initialize Variational Parameters for the Global Rule (theta)
  # theta_mu: [intercept, slope], theta_prec: uncertainty matrix
  theta_mu <- c(1.0, 0.0) 
  theta_prec <- diag(0.1, 2) 
  
  beta_i <- rep(1.0, n)
  X_mat <- cbind(1, covariate_x) # Design matrix for the prior rule
  
  # --- 2. Main Variational Coordinate Ascent Loop (The EM Loop) ---
  for (iter in 1:max_iter) {
    gamma_old <- gamma
    
    # STEP A: Update Gamma (Gene states)
    # Uses current beta_i and neighbor influence
    neighbor_influence <- as.numeric(W_norm %*% gamma)
    eta <- log_lr + alpha + beta_i * neighbor_influence
    gamma <- 1 / (1 + exp(-eta))
    gamma <- pmin(pmax(gamma, 1e-10), 1 - 1e-10) # Numerical stability
    
    # STEP B: Update Beta_i (Local Sensitivity)
    # Propagate neighbor uncertainty using the variance term
    # Var(sum w_ij z_j) = sum w_ij^2 * gamma_j(1-gamma_j)
    E_neighbor_sq <- as.numeric(W_norm_sq %*% (gamma * (1 - gamma)))
    
    # Informed prior mean derived from the current global rule theta
    informed_prior_means <- as.numeric(X_mat %*% theta_mu)
    
    # Denominator (Precision) includes prior variance and network noise
    local_precision <- (1 / prior_var + neighbor_influence^2 + E_neighbor_sq)
    beta_i <- (informed_prior_means / prior_var + gamma * neighbor_influence) / local_precision
    beta_i <- pmax(beta_i, 0) # Physical constraint: sensitivity must be non-negative
    
    # STEP C: Update Theta (Global Prior Rule)
    # Step-by-step Bayesian update: accumulating information about the covariate effect
    # new_prec = prior_prec + (X'X)/prior_var
    theta_prec_new <- theta_prec + (t(X_mat) %*% X_mat) / prior_var
    
    # Update Mean: Precision-weighted combination of current belief and new evidence
    theta_mu <- solve(theta_prec_new) %*% (theta_prec %*% theta_mu + t(X_mat) %*% beta_i / prior_var)
    theta_prec <- theta_prec_new
    
    # STEP D: Update Alpha (Global baseline)
    alpha <- update_alpha_stable(gamma, p_values, mixture, beta_i, W)
    
    # --- 3. Convergence Check ---
    delta <- max(abs(gamma - gamma_old))
    if (delta < tol) break
  }
  
  return(beta_i)
}