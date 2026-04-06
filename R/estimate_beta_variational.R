#' @title Estimate βᵢ using variational inference
#'
#' @description Full Bayesian inference for gene-specific network sensitivity
#' parameters βᵢ with prior regularization. This is the most accurate method
#' but computationally intensive.
#'
#' @param p_values Numeric vector of p-values (length n).
#' @param W Normalized adjacency matrix (n x n).
#' @param mixture List containing mixture components (f0, f1, pi0, pi1) from
#'   \code{estimate_mixture_components}.
#' @param prior_mean Prior mean for βᵢ.
#' @param prior_var Prior variance for βᵢ.
#' @param max_iter Maximum number of EM iterations (default 50).
#' @param tol Convergence tolerance based on change in gamma (default 1e-4).
#' @param verbose If TRUE, print progress messages.
#'
#' @return Numeric vector of estimated βᵢ (length n).
#'
#' @keywords internal
estimate_beta_variational <- function(p_values, W, mixture,
                                      prior_mean, prior_var,
                                      max_iter = 50, tol = 1e-4,
                                      verbose = TRUE) {
  
  n <- length(p_values)
  
  if (verbose) cat("  Running variational inference for βᵢ...\n")
  
  # Initialize
  gamma <- 1 - compute_local_fdr_initial(p_values, mixture)
  alpha <- 0
  beta_i <- rep(prior_mean, n)
  beta_i_var <- rep(prior_var, n)  # Posterior variance
  
  # Precompute for efficiency
  degree <- rowSums(W > 0)
  
  # EM iterations
  for (iter in 1:max_iter) {
    gamma_old <- gamma
    
    # --- E-step: Update gamma with current beta_i ---
    neighbor_influence <- W %*% gamma 
    log_lr <- compute_log_likelihood_ratio(p_values, mixture)
    
    eta <- log_lr + alpha + beta_i * neighbor_influence
    gamma_new <- 1 / (1 + exp(-eta))
    gamma_new <- pmin(pmax(gamma_new, 1e-10), 1 - 1e-10)
    gamma <- gamma_new
    
    # --- M-step: Update beta_i ---
    for (i in 1:n) {
      if (degree[i] > 0) {
        # Compute sufficient statistics
        neighbor_sum <- sum(W[i, ] * gamma)
        #neighbor_sum <- sum(W[i, ] * gamma)
        
        # Bayesian update for beta_i
        # Prior: N(prior_mean, prior_var)
        # Likelihood: gamma[i] ~ Bernoulli(sigmoid(...))
        # Use Laplace approximation for variational update
        
        # Current estimate of linear predictor without beta_i term
        eta_without_beta <- log_lr[i] + alpha
        
        # Compute gradient and Hessian for beta_i
        # Using logistic regression likelihood
        
        # Expected value of neighbor influence squared
        E_neighbor_sq <- sum((W[i, ])^2 * gamma * (1 - gamma))
        
        # Update rule (simplified variational Bayes)
        # beta_i posterior mean
        beta_i[i] <- (prior_mean / prior_var + gamma[i]*neighbor_sum) / 
          (1/prior_var + neighbor_sum^2 + E_neighbor_sq + 1e-10)
        
        # beta_i posterior variance (approximate)
        #beta_i_var[i] <- 1 / (1/prior_var + neighbor_sum^2 + E_neighbor_sq + 1e-10)
        beta_i_var[i] <- 1 / (1/prior_var + neighbor_sum^2 + E_neighbor_sq)
        # Ensure positivity
        beta_i[i] <- max(beta_i[i], 0)
      } else {
        # Isolated genes: shrink toward prior
        beta_i[i] <- prior_mean * 0.1
        beta_i_var[i] <- prior_var
      }
    }
    
    # Update alpha (global baseline)
    alpha <- update_alpha_stable(gamma, p_values, mixture, beta_i, W)
    
    # Check convergence
    delta <- max(abs(gamma - gamma_old))
    if (verbose && iter %% 1 == 0) {
      cat(sprintf("    Variational iter %d: delta = %.2e, mean(βᵢ) = %.3f\n",
                  iter, delta, mean(beta_i)))
    }
    
    if (delta < tol) break
  }
  
  if (verbose) cat(sprintf("Variational inference converged in %d iterations\n",iter))
  
  return(beta_i)
}
