# ============================================================================
# Network Density Selection via Power Law Fit
# ============================================================================

#' @title WGCNA-style power law density selection
#'
#' @description Selects an optimal density threshold for a weighted network
#'   such that the degree distribution follows a power law. This mimics the
#'   scale-free topology criterion used in WGCNA.
#'
#' @param W_raw Raw symmetric weight matrix (n x n).
#' @param min_density Minimum network density to consider (default 0.01).
#' @param max_density Maximum network density to consider (default 0.3).
#' @param n_bins Number of density thresholds to test (default 20).
#' @param Rsquared_cutoff Minimum R² for a fit to be considered good (default 0.85).
#'
#' @return A list containing:
#'   \item{W_sparse}{Sparse binary adjacency matrix after thresholding.}
#'   \item{optimal_density}{Selected density value.}
#'   \item{fit_results}{Data frame with evaluation metrics for each candidate density.}
#'
#' @export
power_law_selection <- function(W_raw,
                                min_density = 0.01,
                                max_density = 0.3,
                                n_bins = 20,
                                Rsquared_cutoff = 0.85) {
  
  n <- nrow(W_raw)
  cat("Power law density selection...\n")
  
  # Candidate densities
  density_candidates <- seq(min_density, max_density, length.out = n_bins)
  fit_results <- data.frame(
    density = density_candidates,
    Rsquared = NA_real_,
    alpha_est = NA_real_,
    n_edges = NA_integer_
  )
  
  # Evaluate each candidate
  for (i in seq_along(density_candidates)) {
    density <- density_candidates[i]
    fit <- evaluate_power_law_fit(W_raw, density)
    fit_results[i, "Rsquared"] <- fit$Rsquared
    fit_results[i, "alpha_est"] <- fit$alpha_est
    fit_results[i, "n_edges"] <- fit$n_edges
  }
  
  # Find optimal density
  optimal_idx <- find_optimal_density(fit_results,
                                      target_alpha = 1.5,           # typical scale-free exponent
                                      Rsquared_cutoff = Rsquared_cutoff)
  optimal_density <- fit_results$density[optimal_idx]
  
  cat(sprintf("Optimal density: %.3f (R²=%.3f, α=%.2f)\n",
              optimal_density,
              fit_results$Rsquared[optimal_idx],
              fit_results$alpha_est[optimal_idx]))
  
  # Create sparse matrix at optimal density
  W_sparse <- apply_density_threshold(W_raw, optimal_density)
  
  return(list(
    W_sparse = W_sparse,
    optimal_density = optimal_density,
    fit_results = fit_results
  ))
}


#' @title Evaluate power law fit for a given density threshold
#'
#' @description Applies a density threshold to a network and fits a power law
#'   to the degree distribution, returning goodness-of-fit metrics.
#'
#' @param W_raw Raw weight matrix.
#' @param density Target network density (proportion of edges retained).
#' @param target_alpha Target power law exponent (default 1.5).
#'
#' @return A list with components:
#'   \item{Rsquared}{Coefficient of determination of the log-log linear fit.}
#'   \item{alpha_est}{Estimated power law exponent (negative slope).}
#'   \item{alpha_diff}{Absolute difference from target_alpha.}
#'   \item{n_edges}{Number of edges after thresholding.}
#'   \item{intercept}{Intercept of the log-log fit.}
#'   \item{rmse}{Root mean squared error of the fit.}
#'   \item{fit_valid}{Logical indicating whether fit meets basic criteria.}
#'   \item{score}{Combined quality score (Rsquared / (1+alpha_diff+rmse)).}
#'   \item{degree_values}{Midpoints of degree bins used for fitting.}
#'   \item{degree_freq}{Frequencies of degree bins.}
#'   \item{log_fit}{lm object of the log-log fit.}
#'   \item{message}{Status message.}
#'
#' @export
evaluate_power_law_fit <- function(W_raw, density, target_alpha = 1.5) {
  n <- nrow(W_raw)
  
  # Base result template
  base_result <- list(
    Rsquared = 0,
    alpha_est = 0,
    alpha_diff = Inf,
    n_edges = 0,
    intercept = 0,
    rmse = Inf,
    fit_valid = FALSE,
    score = 0,
    degree_values = numeric(0),
    degree_freq = numeric(0),
    log_fit = NULL,
    message = "Initialization failed"
  )
  
  # Input checks
  if (n < 10) {
    base_result$message <- "Matrix too small for reliable fitting"
    return(base_result)
  }
  if (density <= 0 || density > 0.5) {
    base_result$message <- "Invalid density value (must be in (0, 0.5])"
    return(base_result)
  }
  
  tryCatch({
    # Apply density threshold
    W_thresholded <- apply_density_threshold(W_raw, density)
    
    # Convert to binary adjacency matrix
    A_binary <- ifelse(W_thresholded > 0, 1, 0)
    n_edges <- sum(A_binary) / 2   # because symmetric, count each edge once
    base_result$n_edges <- n_edges
    
    if (n_edges < n) {
      base_result$message <- "Insufficient edges after thresholding"
      return(base_result)
    }
    
    # Compute degree distribution
    degree_dist <- rowSums(A_binary)
    non_zero_degrees <- degree_dist[degree_dist > 0]
    if (length(non_zero_degrees) < 0.5 * n) {
      base_result$message <- "Too many isolated nodes"
      return(base_result)
    }
    
    # Bin degrees to reduce noise
    deg_hist <- hist(non_zero_degrees, breaks = "FD", plot = FALSE)
    deg_values <- deg_hist$mids[deg_hist$counts > 0]
    deg_freq <- deg_hist$counts[deg_hist$counts > 0]
    
    if (length(deg_values) < 3) {
      base_result$message <- "Insufficient degree bins for fitting"
      return(base_result)
    }
    
    # Log transform
    log_deg <- log(deg_values)
    log_freq <- log(deg_freq)
    
    # Remove infinite values
    valid_idx <- is.finite(log_deg) & is.finite(log_freq)
    if (sum(valid_idx) < 3) {
      base_result$message <- "Invalid values after log transformation"
      return(base_result)
    }
    log_deg <- log_deg[valid_idx]
    log_freq <- log_freq[valid_idx]
    
    # Linear fit: log(freq) ~ log(degree)
    power_fit <- lm(log_freq ~ log_deg)
    
    # Extract metrics
    Rsquared <- summary(power_fit)$r.squared
    alpha_est <- -coef(power_fit)[2]   # negative slope
    intercept <- coef(power_fit)[1]
    residuals <- resid(power_fit)
    rmse <- sqrt(mean(residuals^2))
    
    # Validate fit quality
    fit_valid <- (Rsquared >= 0.7) &&
      (alpha_est > 1.5) && (alpha_est < 4) &&
      (rmse < 1.0)
    
    alpha_diff <- abs(alpha_est - target_alpha)
    score <- Rsquared / (1 + alpha_diff + rmse)
    
    return(list(
      Rsquared = Rsquared,
      alpha_est = alpha_est,
      alpha_diff = alpha_diff,
      n_edges = n_edges,
      intercept = intercept,
      rmse = rmse,
      fit_valid = fit_valid,
      score = score,
      degree_values = deg_values,
      degree_freq = deg_freq,
      log_fit = power_fit,
      message = if (fit_valid) "Good power law fit" else
        paste("Poor fit: R²=", round(Rsquared, 3), ", α=", round(alpha_est, 3))
    ))
    
  }, error = function(e) {
    base_result$message <- paste("Fitting error:", e$message)
    return(base_result)
  })
}


# ============================================================================
# Helper Functions (Internal)
# ============================================================================

#' @title Find optimal density from fit results
#'
#' @param fit_results Data frame with columns: density, Rsquared, alpha_est.
#' @param target_alpha Desired power law exponent (default 1.5).
#' @param Rsquared_cutoff Minimum acceptable R² (default 0.85).
#'
#' @return Index of the optimal density row.
#' @keywords internal
find_optimal_density <- function(fit_results, target_alpha = 1.5, Rsquared_cutoff = 0.85) {
  # First, try to find densities with R² above cutoff
  good_fits <- which(fit_results$Rsquared >= Rsquared_cutoff)
  if (length(good_fits) > 0) {
    # Among those, pick the one with alpha closest to target
    alpha_diff <- abs(fit_results$alpha_est[good_fits] - target_alpha)
    best_idx <- good_fits[which.min(alpha_diff)]
    return(best_idx)
  } else {
    # Otherwise, just pick the highest R²
    return(which.max(fit_results$Rsquared))
  }
}


#' @title Apply density threshold to a weighted matrix
#'
#' @description Keeps the strongest edges (by absolute weight) such that the
#'   resulting network has density approximately equal to the target.
#'
#' @param W_raw Symmetric weight matrix (n x n).
#' @param density Target proportion of edges (excluding diagonal).
#'
#' @return Sparse matrix (same class as W_raw) with thresholded weights.
#' @keywords internal
apply_density_threshold <- function(W_raw, density) {
  n <- nrow(W_raw)
  # Number of possible undirected edges (excluding diagonal)
  total_edges <- n * (n - 1) / 2
  target_edges <- round(density * total_edges)
  
  # Extract upper triangle weights (excluding diagonal)
  upper_tri <- which(upper.tri(W_raw), arr.ind = TRUE)
  weights <- W_raw[upper_tri]
  
  # Sort by absolute weight descending
  order_idx <- order(abs(weights), decreasing = TRUE)
  kept_idx <- order_idx[1:min(target_edges, length(order_idx))]
  
  # Build sparse adjacency matrix
  W_sparse <- matrix(0, n, n)
  for (k in kept_idx) {
    i <- upper_tri[k, 1]
    j <- upper_tri[k, 2]
    w <- weights[k]
    W_sparse[i, j] <- w
    W_sparse[j, i] <- w
  }
  
  return(W_sparse)
}