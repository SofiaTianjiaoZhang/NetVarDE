#' @title Estimate mixture components (null and alternative) from p-values
#'
#' @description Estimates the null (uniform) and alternative distributions from
#'   a vector of p-values using either kernel density estimation (KDE) on the
#'   probit scale or a Beta distribution fit via EM. Also estimates the proportion
#'   of null hypotheses (pi0) using a robust median-based method.
#'
#' @param p_values Numeric vector of p-values (should be in [0,1]).
#' @param method Estimation method: \code{"kde"} (kernel density on probit scale)
#'   or \code{"beta"} (Beta distribution via EM). Default \code{"kde"}.
#' @param verbose If TRUE, print progress messages.
#'
#' @return A list with components:
#'   \item{f0}{Function returning the null density (uniform on [0,1]).}
#'   \item{f1}{Function returning the alternative density (normalized).}
#'   \item{pi0}{Estimated proportion of null p-values.}
#'   \item{pi1}{Estimated proportion of alternative p-values (1 - pi0).}
#'   \item{method}{The method used.}
#'
#' @keywords interna
#' 
#' 
estimate_mixture_components <- function(p_values, method = "kde", verbose = TRUE) {
  
  # Input validation
  if (any(p_values < 0 | p_values > 1)) {
    stop("p-values must be in [0, 1]")
  }
  
  n <- length(p_values)
  
  # 1. Robust pi0 estimation (using multiple lambda)）
  if (verbose) cat("  Estimating pi0...\n")
  
  pi0_robust <- function(p) {
    lambda_seq <- seq(0.05, 0.95, by = 0.05)
    pi0_est <- sapply(lambda_seq, function(l) {
      mean(p > l) / (1 - l)
    })
    # use median for robustness
    median(pi0_est, na.rm = TRUE)
  }
  
  pi0 <- min(0.99, max(0.01, pi0_robust(p_values)))  # clamp to reasonable range
  
  # 2. Null distribution (uniform)
  f0 <- function(p) {
    ifelse(p >= 0 & p <= 1, 1, 0)
  }
  
  # 3. Alternative distribution estimation 
  if (method == "kde") {
    if (verbose) cat("  Estimating f1 via KDE on probit scale...\n")
    
    eps <- 1e-8
    p_adj <- pmin(pmax(p_values, eps), 1 - eps)
    z_scores <- qnorm(p_adj)
    
    # Kernel density estimation on z-scale
    bw <- bw.SJ(z_scores)  # Sheather-Jones bandwidth selector
    kde <- density(z_scores, bw = bw, n = 2048, 
                   from = qnorm(1e-6), to = qnorm(1 - 1e-6))
    
    f_z_mix <- approxfun(kde$x, kde$y, rule = 2)
    
    # Alternative density on z-scale: f1_z = (f_mix - pi0 * phi) / (1 - pi0)
    f1_z <- function(z) {
      mix_dens <- f_z_mix(z)
      null_dens <- dnorm(z)  # 零假设下 z ~ N(0,1)
      alt_dens <- (mix_dens - pi0 * null_dens) / (1 - pi0)
      pmax(0, alt_dens)  # 非负约束
    }
    
    # Transform to p-scale: f1(p) = f1_z(qnorm(p)) * |dz/dp|
    # dz/dp = 1 / dnorm(z) because p = Phi(z) => dp/dz = phi(z)
    f1 <- function(p) {
      p_safe <- pmin(pmax(p, eps), 1 - eps)
      z <- qnorm(p_safe)
      dz_dp <- 1 / dnorm(z)  # Jacobian 项
      f1_val <- f1_z(z) * dz_dp
      
      # 归一化确保在 [0,1] 上积分为 1
      f1_val
    }
    
  } else if (method == "beta") {
    if (verbose) cat("  Estimating f1 via Beta distribution...\n")
    
    # EM algorithm for Beta distribution (with optional zero-weight protection)）
    fit_beta_em <- function(p, pi0, max_iter = 100, tol = 1e-6) {
      alpha <- 0.5
      beta <- 5
      
      for (iter in 1:max_iter) {
        # E-step: posterior probability of being alternative
        f0_val <- 1  
        f1_val <- dbeta(p, alpha, beta)
        post_prob <- (1 - pi0) * f1_val / (pi0 * f0_val + (1 - pi0) * f1_val)
        
        # M-step: weighted method of moments
        w <- post_prob
        m <- sum(w * p) / sum(w)
        v <- sum(w * (p - m)^2) / sum(w)
     
        if (v < m * (1 - m)) {
          term <- m * (1 - m) / v - 1
          alpha_new <- max(0.1, m * term)
          beta_new <- max(0.1, (1 - m) * term)
        } else {
          alpha_new <- alpha
          beta_new <- beta
        }
        
        # Check convergence
        if (abs(alpha_new - alpha) < tol && abs(beta_new - beta) < tol) {
          alpha <- alpha_new
          beta <- beta_new
          break
        }
        
        alpha <- alpha_new
        beta <- beta_new
      }
      
      return(list(alpha = alpha, beta = beta))
    }
    
    params <- fit_beta_em(p_values, pi0)
    f1 <- function(p) dbeta(p, params$alpha, params$beta)
  }
  
  # 4. Numerical normalization of f1 (to ensure it integrates to 1)
  if (verbose) cat("  Normalizing f1...\n")
  
  integrate_f1 <- function() {
    # Focus on [0, 0.1] where alternative density is typically concentrated
    pts1 <- seq(1e-6, 0.1, length.out = 500)
    pts2 <- seq(0.1, 1 - 1e-6, length.out = 200)
    test_pts <- c(pts1, pts2)
    
    f1_vals <- f1(test_pts)
    
    # Clip extreme values to avoid integration issues
    f1_vals <- pmin(f1_vals, quantile(f1_vals, 0.99, na.rm = TRUE))
    f1_vals[!is.finite(f1_vals)] <- 0
    
    # Simple trapezoidal integration
    n <- length(test_pts)
    h <- diff(test_pts)
    integral <- sum((f1_vals[-n] + f1_vals[-1]) * h / 2)
    
    integral
  }
  
  integral <- integrate_f1()
  
  if (integral > 0) {
    f1_final <- function(p) f1(p) / integral
  } else {
    warning("f1 integrates to 0, using default Beta(0.5, 5)")
    f1_final <- function(p) dbeta(p, 0.5, 5)
  }
  
  return(list(
    f0 = f0,
    f1 = f1_final,
    pi0 = pi0,
    pi1 = 1 - pi0,
    method = method
  ))
}
