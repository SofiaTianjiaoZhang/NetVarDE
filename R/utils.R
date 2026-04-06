# ============================================================================
# Utility Functions for NELFDR
# ============================================================================

#' @title Compute initial local false discovery rate (lfdr)
#'
#' @description Calculates the local false discovery rate from mixture model
#'   estimates. This is used to initialize the posterior probabilities.
#'
#' @param p_values Numeric vector of p-values.
#' @param mixture List with components f0, f1, pi0, pi1 from
#'   \code{estimate_mixture_components}.
#'
#' @return Numeric vector of lfdr values (same length as p_values).
#'
#' @keywords internal
compute_local_fdr_initial <- function(p_values, mixture) {
  f0_vals <- mixture$f0(p_values)
  f1_vals <- mixture$f1(p_values)
  lfdr <- mixture$pi0 * f0_vals / (mixture$pi0 * f0_vals + mixture$pi1 * f1_vals)
  lfdr[is.na(lfdr)] <- 1
  lfdr <- pmin(pmax(lfdr, 0), 1)
  return(lfdr)
}


#' @title Compute log likelihood ratio log(f1/f0)
#'
#' @description Computes the log-likelihood ratio for each p-value under the
#'   alternative versus null distributions. Handles zeros and infinite values
#'   by clamping to a reasonable range.
#'
#' @param p_values Numeric vector of p-values.
#' @param mixture List with components f0, f1 from \code{estimate_mixture_components}.
#'
#' @return Numeric vector of log-likelihood ratios.
#'
#' @keywords internal
compute_log_likelihood_ratio <- function(p_values, mixture) {
  f0_vals <- mixture$f0(p_values)
  f1_vals <- mixture$f1(p_values)
  f0_vals[f0_vals == 0] <- .Machine$double.eps
  f1_vals[f1_vals == 0] <- .Machine$double.eps
  log_lr <- log(f1_vals) - log(f0_vals)
  log_lr[is.na(log_lr)] <- 0
  log_lr[is.infinite(log_lr)] <- sign(log_lr[is.infinite(log_lr)]) * 10
  return(log_lr)
}


#' @title Normalize network adjacency matrix
#'
#' @description Applies various normalization methods to a square adjacency
#'   matrix. Ensures symmetry, removes self-loops, and optionally row-normalizes,
#'   symmetrically normalizes, binarizes, or leaves unchanged.
#'
#' @param W Square adjacency matrix (potentially asymmetric, with self-loops).
#' @param method Normalization method: 
#'   \itemize{
#'     \item \code{"row"}: row-normalize so each row sums to 1.
#'     \item \code{"sym"}: symmetric normalization \eqn{D^{-1/2} W D^{-1/2}}.
#'     \item \code{"none"}: no normalization (only symmetrization and diag removal).
#'     \item \code{"binary"}: convert to binary (presence/absence) matrix.
#'   }
#' @param verbose If TRUE, print warnings about isolated nodes.
#'
#' @return Normalized matrix of same dimensions as W.
#'
#' @keywords internal
normalize_network_matrix <- function(W, method = "row", verbose = TRUE) {
  # 1. Symmetrize if not symmetric
  if (!isSymmetric(W)) {
    if (verbose) cat("  Making matrix symmetric...\n")
    W <- (W + t(W)) / 2
  }
  
  # 2. Remove self-loops
  diag(W) <- 0
  
  # 3. Normalize according to method
  if (method == "row") {
    row_sums <- rowSums(W)
    zero_rows <- row_sums == 0
    if (any(zero_rows)) {
      if (verbose) cat(sprintf("  Warning: %d isolated nodes (zero row sum)\n", sum(zero_rows)))
      row_sums[zero_rows] <- 1  # avoid division by zero
    }
    W_norm <- W / row_sums
    return(W_norm)
    
  } else if (method == "sym") {
    D <- diag(1 / sqrt(rowSums(W) + 1e-10))
    W_norm <- D %*% W %*% D
    diag(W_norm) <- 0
    return(W_norm)
    
  } else if (method == "none") {
    return(W)
    
  } else if (method == "binary") {
    W_bin <- (W > 0) * 1
    diag(W_bin) <- 0
    return(W_bin)
    
  } else {
    stop("Unknown normalization method. Choose from: 'row', 'sym', 'none', 'binary'")
  }
}