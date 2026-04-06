#' @title Update the global intercept α with fixed βᵢ (numerically stable)
#'
#' @description Performs a Newton-Raphson update for the global intercept α
#' given the current posterior probabilities γ and gene-specific βᵢ. The update
#' uses numerically stable sigmoid approximations to avoid overflow.
#'
#' @param gamma Numeric vector of current posterior probabilities (length n).
#' @param p_values Numeric vector of p-values (length n).
#' @param mixture List with components f0, f1, pi0, pi1 from
#'   \code{estimate_mixture_components}.
#' @param beta_i Numeric vector of gene-specific βᵢ (length n).
#' @param W Normalized adjacency matrix (n x n).
#' @param max_iter Maximum number of Newton iterations (default 20).
#' @param tol Convergence tolerance based on change in α (default 1e-6).
#'
#' @return Updated α value (numeric scalar).
#'
#' @keywords internal
update_alpha_stable <- function(gamma, p_values, mixture, beta_i, W,
                                max_iter = 20, tol = 1e-6) {
  
  n <- length(gamma)
  
  # Compute neighbor influence once per call
  
  neighbor_influence <- W %*% gamma
  
  # Precompute log-likelihood ratio
  log_lr <- compute_log_likelihood_ratio(p_values, mixture)
  
  # Initialize
  alpha_current <- 0
  alpha_prev <- NA
  
  for (iter in 1:max_iter) {
    # Current linear predictor
    eta <- log_lr + alpha_current + beta_i * neighbor_influence
    
    # Numerically stable sigmoid (piecewise)
    mu <- numeric(n)
    large_pos <- eta > 10
    large_neg <- eta < -10
    moderate <- !(large_pos | large_neg)
    
    mu[moderate] <- 1 / (1 + exp(-eta[moderate]))
    
    mu[large_pos] <- 1 - exp(-eta[large_pos])
    mu[large_pos] <- pmin(pmax(mu[large_pos], 1e-10), 1 - 1e-10)
    
    mu[large_neg] <- exp(eta[large_neg])
    mu[large_neg] <- pmin(pmax(mu[large_neg], 1e-10), 1 - 1e-10)
    
    # Gradient and Hessian
    gradient <- sum(gamma - mu)
    hessian <- -sum(mu * (1 - mu))
    
    # Regularization to avoid division by zero
    if (abs(hessian) < 1e-10) {
      hessian <- -1e-10 * sign(hessian)
    }
    
    # 牛顿更新，带步长缩减
    step <- -gradient / hessian
    
    # Newton step
    max_step <- 2.0
    if (abs(step) > max_step) {
      step <- sign(step) * max_step
    }
    
    alpha_new <- alpha_current + step
    
    # Limit step size for stability
    if (!is.na(alpha_prev) && abs(alpha_new - alpha_prev) < tol) {
      return(alpha_new)
    }
    
    alpha_prev <- alpha_current
    alpha_current <- alpha_new
  }
  
  return(alpha_current)
}