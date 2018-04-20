#' @title Local inidcator of spatial association
#' @description \code{lisa} is a function to estimate the local indicators of spatial association. The function assumes univariate data at each location. For multivariate data use \code{\link{lisa.nc}}
#' @param x vector of length n representing the x coordinates (or latitude; see latlon).
#' @param y vector of length n representing the y coordinates (or longitude).
#' @param z vector of n representing the observation at each location.
#' @param neigh neighborhood size.
#' @param resamp number of resamples under the NULL to generate p-values
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "lisa" is returned, consisting of the following components: 
#' \item{correlation}{the autocorrelation within the neighborhood (neigh) of each observation measured using Moran's I.}
#' \item{p}{the permutation two-sided p-value for each observation.}
#' \item{mean}{the mean of the observations inside each neighborhooddistance within each neighborhood.}
#' \item{n}{the number of observations within each neighborhood.}
#' \item{dmean}{the actual mean distance within each neighborhood.}
#' \item{z}{the original observations}
#' \item{coord}{a list with the x and y coordinates.}
#' @details This is the function to estimate the local indicators of spatial association modified form Anselin (1995). The statistic is the average autocorrelation within a neighborhood.
#' @references Anselin, L. 1995. Local indicators of spatial association - LISA. Geographical Analysis 27:93-115. \url{https://doi.org/10.1111/j.1538-4632.1995.tb00338.x}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{plot.lisa}} 
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[, 2]
#' 
#' # z data from an exponential random field
#' z <- rmvn.spa(x = x, y = y, p = 2, method = "gaus")
#' 
#' # lisa analysis
#' fit1 <- lisa(x = x, y = y, z = z, neigh = 3, resamp = 500)
#' \dontrun{plot(fit1, negh.mean=FALSE)}
#' @keywords spatial
#' @export
################################################################################
lisa <- function(x, y, z, neigh, resamp = 1000, latlon = FALSE, quiet = FALSE) {
  ##############################################################################
  # lisa is a function to estimate the local indicators of spatial association.
  ##############################################################################
  if (!is.null(dim(z))) {
    stop("\n z is multivariate. Use lisa.nc()")
  }
  
  n <- length(z)
  # then generating geographic distances
  if (latlon) {
    # these are geographic distances from lat-lon coordinates
    dmat <- gcdist(x, y)
  } else{
    dmat <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  }
  
  zscal <- (scale(z, center = TRUE, scale = TRUE)[, 1])/(sqrt((n - 1)/n))
  zx <- t(outer(zscal, zscal))
  dkl <- ifelse(dmat > 0 & dmat < neigh, 1, NA)	# flaggs obs within the neigh
  
  # calculate mean value of z within neigh
  negz <- dkl*z
  diag(negz) <- z
  zmean <- apply(negz, 2, mean, na.rm = TRUE)
  
  nlok <- apply(dkl, 2, sum, na.rm = TRUE)
  dmean <- apply(dmat*dkl, 2, mean, na.rm = TRUE)
  moran <- apply(zx*dkl, 2, mean, na.rm = TRUE)
  
  p <- NULL
  if (resamp > 0) {
    perm <- matrix(NA, nrow = resamp, ncol = n)
    
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      trekk <- sample(1:n)
      zx2 <- zx[trekk, trekk]
      perm[i, ] <- apply(zx2*dkl, 2, mean, na.rm = TRUE)
    }
    p <- (apply(moran < t(perm), 1, sum))/(resamp + 1)
    p <- apply(cbind(p, 1 - p), 1, min) + 1/(resamp + 1)
  }
  res <- list(correlation = moran, p = p, mean = zmean, dmean = dmean, n = nlok, 
              z = z, coord = list(x = x, y = y), call = deparse(match.call()))
  
  class(res) <- "lisa"
  res
}

#' @title Plots local indicators of spatial association
#' @description `plot' method for class "lisa".
#' @param x an object of class "lisa", ususally, as a result of a call to \code{\link{lisa}}.
#' @param neigh.mean If TRUE, size of symbols represents average observation in each neighborhood; If FALSE, size of symbols represents the original observation
#' @param add If TRUE, a lisa-plot will be added to a pre-existing plot.
#' @param inches scales the size of the symbols
#' @param \dots other arguments
#' @return A bubble-plot of observations against spatial coordinates is produced. Below mean values are signified by squares. Above mean values are signified by squares. 
#' 
#'   If a permutation test was performed, observations for which the associated LISA statistic is positive and significant at a nominal (two-sided) 5\%-level will be respresented by filled symbols and non-significant values by open symbols. Thus spatial hot-spots are represented by red filled circles and cold-spots by black filled squares.
#' @seealso \code{\link{lisa}}, \code{\link{lisa.nc}}
#' @keywords spatial
#' @export
################################################################################
plot.lisa <- function(x, neigh.mean = FALSE, add = FALSE, inches = 0.2, ...) {
  ##############################################################################
  xx <- x
  if (neigh.mean) {
    z <- xx$mean
  } else {
    z <- xx$z
  }
  
  x <- xx$coord$x
  y <- xx$coord$y
  
  if (add == FALSE) {
    plot(x, y, type = "n", ...)
  }
  sel <- is.finite(z)
  x <- split(x, z - mean(z, na.rm = TRUE) > 0)
  y <- split(y, z - mean(z, na.rm = TRUE) > 0)
  sel <- split(sel, z - mean(z, na.rm = TRUE) > 0)
  z2 <- split(z - mean(z, na.rm = TRUE), z - mean(z, na.rm = TRUE) > 0)
  
  bgc <- rep(0, length(z))
  bgc <- split(bgc, (z - mean(z, na.rm = TRUE)) > 0)
  
  if (!is.null(xx$p)) {
    bgc <- rep(0, length(z))
    bgc <- ifelse(xx$p < 0.025, 2, 0)
    bgc[xx$p < 0.025 & (z - mean(z, na.rm = TRUE)) > 0] <- 1
    bgc <- split(bgc, (z - mean(z, na.rm = TRUE)) > 0)
  }
  
  if (!is.null(length(z2[[1]][sel[[1]]]))) {
    symbols(x[[1]][sel[[1]]], y[[1]][sel[[1]]], circles = -z2[[1]][sel[[1]]], 
            inches = inches, add = TRUE, fg = 2, bg = bgc[[1]][sel[[1]]])
  }
  
  if (!is.null(length(z2[[1]][sel[[2]]]))) {
    symbols(x[[2]][sel[[2]]], y[[2]][sel[[2]]], squares = z2[[2]][sel[[2]]], 
            inches = inches, add = TRUE, fg = 1, bg = bgc[[2]][sel[[2]]])
  }
}

#' @title Non-centered inidcators of spatial association
#' @description \code{lisa.nc} is a function to estimate the (noncentred) multivariate local indicators of spatial association. The function requires multiple observations at each location. For single observations at each location use \code{lisa}.
#' @param x vector of length n representing the x coordinates (or latitude; see latlon).
#' @param y vector of length n representing the y coordinates (or longitude).
#' @param z a matrix of dimension n x p representing p (>1) observation at each location.
#' @param neigh neighborhood size.
#' @param resamp number of resamples under the NULL to generate p-values
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "lisa" is returned, consisting of the following components:
#' \item{correlation}{the mean correlation within the neighborhood (neigh).}
#' \item{p}{the permutation two-sided p-value for each distance-class.}
#' \item{n}{the number of pairs within each neighborhood.}
#' \item{dmean}{the actual mean of distance within each neighborhood.}
#' \item{coord}{a list with the x and y coordinates.}
#' @details This is the function to estimate the (non-centered) local indicators of spatial 
#' association modified form Anselin (1995). 'correlation' is the average correlation within 
#' a neighborhood. The function requires multiple observations at each location.
#' 
#'   Missing values are allowed -- values are assumed missing at random, and pairwise complete observations will be used.
#' @references Anselin, L. 1995. Local indicators of spatial association - LISA. Geographical Analysis 27:93-115. \url{https://doi.org/10.1111/j.1538-4632.1995.tb00338.x}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{lisa}}}
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[,2]
#' 
#' # z data from an exponential random field
#' z <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp")
#'   )
#' 
#' # lisa.nc analysis
#' fit1 <- lisa.nc(x = x, y = y, z = z, neigh = 3)
#' \dontrun{plot(fit1)}
#' @keywords spatial
#' @export
################################################################################
lisa.nc <- function(x, y, z, neigh, na.rm = FALSE, resamp = 1000, latlon = FALSE, 
                  quiet = FALSE) {
  ##############################################################################
  if (is.null(dim(z))) {
    stop("\n z is univariate. Use lisa()")
  }
  
  NAO <- FALSE
  # check for missing values
  if (any(!is.finite(unlist(z)))) {
    if (na.rm) {
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    } else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  n <- dim(z)[1]
  p <- dim(z)[2]
  
  # then generating geographic distances
  if (latlon) {
    # these are geographic distances from lat-lon coordinates
    dmat <- gcdist(x,y)
  } else {
    dmat <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  }
  
  z <- as.matrix(z) + 0
  zx <- cor2(t(z), circ = FALSE)
  
  dkl <- ifelse(dmat > 0 & dmat < neigh, 1, NA)	# flaggs obs within the neigh
  nlok <- apply(dkl, 2, sum, na.rm = TRUE)
  dmean <- apply(dmat*dkl, 2, mean, na.rm = TRUE)
  moran <- apply(zx*dkl, 2, mean, na.rm = TRUE)
  
  p <- NULL
  if (resamp > 0) {
    perm <- matrix(NA, nrow = resamp, ncol = n)
    
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      trekk <- sample(1:n)
      zx2 <- zx[trekk, trekk]
      perm[i, ] <- apply(zx2*dkl, 2, mean, na.rm = TRUE)
    }
    p <- (apply(moran < t(perm), 1, sum))/(resamp + 1)
    p <- apply(cbind(p, 1 - p), 1, min) + 1/(resamp + 1)
  }
  
  res <- list(correlation = moran, p = p, n = nlok, dmean = dmean, 
              coord = list(x = x, y = y), call = deparse(match.call()))
  class(res) <- "lisa"
  res
}
