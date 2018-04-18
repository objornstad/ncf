#' @title Mantel Test
#' @description A simple function to do a permutation-based Mantel test. The data can either be two distance/similarity matrices or (x, y, z) data.
#' @param M1 similarity/distance matrix 1
#' @param M2 similarity/distance matrix 2
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z matrix of dimension n x p representing p observation at each location.
#' @param resamp the number of resamples for the null distribution.
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "Mantel" is returned, consisting of a list with two components:
#' \item{correlation}{the value for the Mantel correlation.}
#' \item{p}{the randomization-based two-sided p-value.}
#' @details Typical usages are
#' \preformatted{
#' mantel.test(M1, M2, x = NULL, y = NULL, z = NULL, resamp = 1000, 
#'             latlon = FALSE, quiet = FALSE)
#' 
#' mantel.test(x, y, z, M1 = NULL, M2 = NULL, resamp = 1000, latlon = FALSE, 
#'             quiet = FALSE)
#' }
#' 
#'   Missing values are treated through pairwise deletion.
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[, 2]
#' # z data from an exponential random field
#' z <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp")
#'   )
#' 
#' # the Mantel test
#' mantel.test(x = x, y = y, z = z[, 1], resamp = 500)
#' @keywords spatial
#' @export
################################################################################
mantel.test <- function(M1 = NULL, M2 = NULL, x = NULL, y = NULL, z = NULL, 
                        resamp = 1000, latlon = FALSE, quiet = FALSE) {
  ##############################################################################
  # mantel.test is a function to calculate the mantel test for two matrices,
  # or for {x, y, z} data.
  ##############################################################################
  if (is.null(M1) & is.null(x)) {
    stop("you must provide either distance/similarity\nmatrices OR vectors of x-/y-coordinates and observations")
  }
  
  if (!is.null(x)) {
    multivar <- !is.null(dim(z)) # test whether z is univariate or multivariate
    
    if (multivar == TRUE) {
      n <- length(z[, 1])
      p <- length(z[1, ])
      z <- as.matrix(z) + 0
      
      M2 <- cor2(t(z), circ = FALSE)
    } else {
      n <- length(z)
      z <- as.vector(z) + 0
      zscal <- (scale(z, center = TRUE, scale = TRUE)[, 1])/(sqrt((n - 1)/n))
      
      M2 <- t(outer(zscal, zscal))
    }
    
    # then generating geographic distances
    if (latlon) {
      # these are geographic distances from lat-lon coordinates
      M1 <- gcdist(x, y)
    } else {
      # these are geographic distances from euclidian coordinates
      M1 <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
    }
  } else {	# if x is null
    n <- dim(M1)[1]
  }
  
  if (resamp != 0) {
    M12 <- M1
  }
  
  M1 <- M1[row(M1) != col(M1)]
  M2 <- M2[row(M2) != col(M2)]
  
  MantelR <- cor2(M1, M2, circ = FALSE)
  
  if (resamp != 0) {
    perm <- rep(NA, resamp)		
    
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (!quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      trekk <- sample(1:n)
      d <- M12[trekk, trekk]
      m <- M2
      d <- d[row(d) != col(d)]
      
      perm[i] <- cor2(d, M2, circ = FALSE)
    }
    p <- (sum(MantelR >= perm))/(resamp + 1)
    p <- min(c(p, 1 - p)) + 1/(resamp + 1)
  } else{ # if resamp == 0
    p <- NA
  }
  
  res <- list(correlation = MantelR, p = p, call = deparse(match.call()))
  
  class(res) <- "Mantel"
  res
}

#' @title Partial Mantel test
#' @description A simple function to calculate permutation-based partial mantel tests for three matrices, the partial mantel test is calculated to test for relationships between M1 and M2 (M3) cotrolling for M3 (M2). syntax and logic follows Legendre and Legendre (1998) pp 557-558.
#' @param M1 similarity/distance matrix 1
#' @param M2 similarity/distance matrix 2
#' @param M3 similarity/distance matrix 3
#' @param resamp the number of resamples for the null distribution.
#' @param method the method to be used for calculating the correlations.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "partial.Mantel" is returned, consisting of a list with two components:
#' \item{MantelR}{the vector of observed Mantel and partial Mantel correlations.}
#' \item{p}{the vector of two-sided p-value under randomization (of M1).}
#' @details Missing values are treated through pairwise deletion. 
#' 
#'   The method must be one of pearson (default), spearman or kendall.
#' @references Legendre, P., and L. Legendre. 1998. Numerical Ecology, 2nd edition. Elsevier, Amsterdam
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{mantel.test}}
#' @examples 
#' # first generate some sample data and dissimilarity matrices
#' x <- rnorm(10)
#' y <- rnorm(10)
#' z <- rnorm(10)
#' M1 <- sqrt(outer(x, x, "-")^2)
#' M2 <- sqrt(outer(y, y, "-")^2)
#' M3 <- sqrt(outer(z, z, "-")^2)
#' 
#' partial.mantel.test(M1 = M1, M2 = M2, M3 = M3, resamp = 500)
#' @keywords spatial
#' @export
################################################################################
partial.mantel.test <- function(M1, M2, M3, resamp = 1000, method = 'pearson', 
                                quiet = FALSE) {
  ##############################################################################
  # partial.mantel.tets is a simple function to calculate Mantel and partial mantel tests for three matrices,
  # the partial mantel test is calculated to test for relationships between M1 and M2 (or M3) cotrolling for M3 (M2).
  # syntax and logic follows Legendre & Legendre (1998) pp 557-558
  ##############################################################################
  # check for missing values
  if (any(!is.finite(M1)) || any(!is.finite(M2)) || any(!is.finite(M3))) {
    warning("Missing values exist; Pairwise deletion will be used")
  }
  
  if (any(outer(c(dim(M1), dim(M2), dim(M3)), c(dim(M1), dim(M2), dim(M3)), "-") != 0)) {
    stop("All matrices must be square and of the same size")
  }
  
  
  n <- dim(M1)[1]
  
  r12 <- cor(M1[upper.tri(M1)], M2[upper.tri(M2)], use = "pairwise.complete.obs", 
             method = method)
  r13 <- cor(M1[upper.tri(M1)], M3[upper.tri(M3)], use = "pairwise.complete.obs", 
             method = method)
  r23 <- cor(M2[upper.tri(M2)], M3[upper.tri(M3)], use = "pairwise.complete.obs", 
             method = method)
  
  r12.3 <- (r12 - r13*r23)/(sqrt(1 - r13^2)*sqrt(1 - r23^2))
  r13.2 <- (r13 - r12*r23)/(sqrt(1 - r12^2)*sqrt(1 - r23^2))
  
  perm <- matrix(NA, ncol = 5, nrow = resamp)
  
  for (i in 1:resamp) {
    whn <- pretty(c(1, resamp), n = 10)
    if (!quiet & any(i == whn)) {
      cat(i, " of ", resamp, "\r")
      flush.console()
    }
    
    trekk <- sample(1:n)
    M1r <- M1[trekk, trekk]
    r12r <- cor(M1r[upper.tri(M1)], M2[upper.tri(M2)], use = "pairwise.complete.obs", 
                method = method)
    r13r <- cor(M1r[upper.tri(M1)], M3[upper.tri(M3)], use = "pairwise.complete.obs", 
                method = method)
    
    perm[i, 1] <- r12r
    perm[i, 2] <- r13r
    
    perm[i, 4] <- (r12r - r13r*r23)/(sqrt(1 - r13r^2)*sqrt(1 - r23^2))
    perm[i, 5] <- (r13r - r12r*r23)/(sqrt(1 - r12r^2)*sqrt(1 - r23^2))
    
    M2r <- M2[trekk, trekk]
    r23r <- cor(M2r[upper.tri(M2)], M3[upper.tri(M3)], use = "pairwise.complete.obs", 
                method = method)
    perm[i, 3] <- r23r
  }
  
  res <- c(r12, r13, r23, r12.3, r13.2)
  names(res) <- c("r12", "r13", "r23", "r12.3", "r13.2")
  
  p <- (apply(t(res >= t(perm)), 2, sum)) / (resamp + 1)
  p <- apply(cbind(p, 1 - p), 1, min) + 1/(resamp + 1)
  
  out <- list(MantelR = res, p = p, call = deparse(match.call()))
  class(out) <- "partial.Mantel"
  
  return(out)
}
