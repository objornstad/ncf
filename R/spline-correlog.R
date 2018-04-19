#' @title Uni- and multivariate spline correlograms
#' @description \code{spline.correlog} is the function to estimate the spline (cross-)correlogram from spatial data. Either univariate or multivariate (time seres) for each site can be used.
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z vector of length n or matrix of dimension n x p representing p observation at each location.
#' @param w an optional second variable with idenitical dimension to z (to estimate cross-correlograms).
#' @param df degrees of freedom for the spline. Default is sqrt(n).
#' @param type takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator.
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param npoints the number of points at which to save the value for the spline function (and confidence envelope / null distribution).
#' @param save If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).
#' @param filter If TRUE, the Fourier filter method of Hall and coworkers is applied to ensure positive semidefiniteness of the estimator.
#' @param fw If filter is TRUE, it may be useful to truncate the function at some distance w sets the truncation distance. When set to zero, no truncation is done.
#' @param max.it the maximum iteration for the Newton method used to estimate the intercepts.
#' @param xmax If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "spline.correlog" is returned, consisting of the following components: 
#' \item{real}{the list of estimates from the data.}
#' \item{$x.intercept}{the lowest value at which the function is = 0. If correlation is initially negative, the distance is given as negative.}
#' \item{$e.intercept}{the lowest value at which the function 1/e.}
#' \item{$y.intercept}{the extrapolated value at x=0 (nugget).}
#' \item{$predicted$x}{the x-axes for the fitted covariance function.}
#' \item{$predcited$y}{the values for the covariance function.}
#' \item{boot}{a list with the analogous output from the bootstrap or null distribution.} 
#' \item{$summary}{gives the full vector of output for the x.intercept, y.intercept, e.intercept, and a quantile summary for the resampling distribution.}
#' \item{$boot}{If save=TRUE, the full raw matrices from the resampling is saved.}
#' \item{max.distance}{the maximum spatial distance considered.}
#' @details If observations are univariate the spline (cross-)correlogram represents the generalization of the spatial (cross-)correlogram; if observations are multivariate the spline (cross-)correlogram represents the generalization of the Mantel (cross-)correlogram.
#' 
#'   The spline (cross-)correlogram differes from the spatial correlogram (and Mantel correlogram) in that it estimated spatial dependence as a continous functions of distance (rather than binning into distance classes). The spline correlogram differs from the nonparametric (cross-)correlation function in that the zero-correlation reference line in the former corresponds to the regionwide correlation reference line in the latter. The x-intercept in the spline correlogram is the distance at which object are no more similar than that expected by-chance-alone across the region. 
#'   
#'   Missing values are allowed -- values are assumed missing at random.
#' @references Bjornstad, O.N. & Falck, W. (2001) Nonparametric spatial covariance functions: estimation and testing. Environmental and Ecological Statistics, 8:53-70. \url{https://doi.org/10.1023/A:1009601932481}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{summary.spline.correlog}}, \code{\link{plot.spline.correlog}}, \code{\link{Sncf}}, \code{\link{spline.correlog2D}}, \code{\link{correlog}}
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[, 2]
#' 
#' # z data from an exponential random field
#' z <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp")
#'   )
#' 
#' # w data from a gaussian random field
#' w <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "gaus"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "gaus")
#'   )
#' 
#' # univariate spline correlogram
#' fit1 <- spline.correlog(x = x, y = y, z = z[, 1], resamp = 100)
#' \dontrun{plot.spline.correlog(fit1)}
#' summary(fit1)
#' 
#' # multivariate spline correlogram
#' fit2 <- spline.correlog(x = x, y = y, z = z, resamp = 100)
#' \dontrun{plot.spline.correlog(fit2)}
#' summary(fit2)
#' 
#' # multivariate spline cross-correlogram
#' fit3 <- spline.correlog(x = x, y = y, z = z, w = w, resamp = 100)
#' \dontrun{plot.spline.correlog(fit3)}
#' summary(fit3)
#' @keywords smooth spatial
#' @export
################################################################################
spline.correlog <- function(x, y, z, w = NULL, df = NULL, type = "boot", resamp = 1000, 
                            npoints = 300, save = FALSE, filter = FALSE, fw = 0, 
                            max.it = 25, xmax = FALSE, latlon = FALSE, na.rm = FALSE, 
                            quiet = FALSE) {
  ##############################################################################
  # spline.correlog is the function to estimate the spline correlogram discussed in 
  # Bjornstad & Falck (2001, Evironmental and Ecological Statistics 8:53-70)
  # Bjornstad et al. (1999, Trends in Ecology and Evolution 14:427-431)
  # Bjornstad et al. (1999, Ecology 80:622-637)
  ##############################################################################
  multivar <- !is.null(dim(z)) # test whether z is univariate or multivariate
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
  
  # This generates the distance matrices
  # first generating the moran distances
  zx <- NULL
  if (multivar == TRUE) {
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z) + 0
    zscal <- (t(apply(z, 2, scale, center = TRUE, scale = TRUE)))/(sqrt((n - 1)/n))
    
    if (is.null(w)) {
      moran <- cor2(t(z), circ = FALSE)
      moran <- moran - mean(moran[lower.tri(moran)], na.rm = TRUE)
    } else {
      w <- as.matrix(w) + 0
      moran <- cor2(t(z), t(w), circ = FALSE)
      moran <- moran - mean(moran, na.rm = TRUE)
    }
  } else {
    n <- length(z)
    z <- as.vector(z) + 0
    zscal <- (scale(z, center = TRUE, scale = TRUE)[, 1])/(sqrt((n - 1)/n))
    
    if (is.null(w)) {
      moran <- t(outer(zscal, zscal))
    } else {
      wscal <- (scale(w, center = TRUE, scale = TRUE)[, 1])/(sqrt((n - 1)/n))
      zw <- c(zscal, wscal)
      moran <- t(outer(zw, zw))[1:n, (n + 1):(2*n)]
    }
  }
  
  if (is.null(df)) {
    df <- sqrt(n)
  }
  
  # then generating geographic distances
  if (latlon) {
    # these are geographic distances from lat-lon coordinates
    xdist <- gcdist(x, y)
  } else {
    # these are geographic distances from euclidian coordinates
    xdist <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  }
  maxdist <- ifelse(!xmax, max(na.omit(xdist)), xmax)
  
  # The spline function
  if (is.null(w)) {
    triang <- lower.tri(xdist)
  } else {
    triang <- is.finite(xdist)
  }
  
  u <- xdist[triang]
  v <- moran[triang]
  sel <- is.finite(u) & is.finite(v)
  u <- u[sel]
  v <- v[sel]
  
  v <- v[u <= maxdist]	
  u <- u[u <= maxdist]
  
  xpoints <- seq(0, maxdist, length = npoints)
  out <- gather(u = u, v = v, w = w, moran = moran, df = df, xpoints = xpoints, 
                filter = filter, fw = fw)
  real <- list(x.intercept = out$xint, e.intercept = out$eint, y.intercept = out$yint,
               predicted = list(x = matrix(out$x, nrow = 1),
                                y = matrix(out$y, nrow = 1)))
  
  # end of spline fit
  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if (resamp != 0) {
    # here is the bootstrapping/randomization
    boot$boot.summary$x.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), 
                      y = matrix(NA, nrow = resamp, ncol = npoints))
    type <- charmatch(type, c("boot", "perm"), nomatch = NA)
    if (is.na(type))
      stop("method should be \"boot\", or \"perm\"")
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (quiet & any(i == whn)) {
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      if (type == 1) {
        trekkx <- sample(1:n, replace = TRUE)
        trekky <- trekkx
      }
      #if(type == 2) {
      #	trekky <- sample(1:n, replace = FALSE)
      #	trekkx <- 1:n
      #}
      xdistb <- xdist[trekkx, trekkx]
      
      if (is.null(w)) {
        triang <- lower.tri(xdist)
      } else {
        triang <- is.finite(xdist)
      }
      
      xdistb <- xdistb[triang]
      moranb <- moran[trekky, trekky][triang]
      if (type == 1 & is.null(w)) {
        moranb <- moranb[!(xdistb == 0)]
        xdistb <- xdistb[!(xdistb == 0)]
      }
      u <- xdistb
      v <- moranb
      sel <- is.finite(u) & is.finite(v)
      u <- u[sel]
      v <- v[sel]
      
      v <- v[u <= maxdist]	
      u <- u[u <= maxdist]
      xpoints <- seq(0, maxdist, length = npoints)
      out <- gather(u = u, v = v, w = w, moran = moranb, df = df, xpoints = xpoints, 
                    filter = filter, fw = fw)
      
      boot$boot.summary$y.intercept[i, 1] <- out$yint
      boot$boot.summary$x.intercept[i, 1] <- out$xint
      boot$boot.summary$e.intercept[i, 1] <- out$eint
      predicted$x[1, ] <- out$x
      predicted$y[i, ] <- out$y
      
    }
    # end of bootstrap loop!
    
    if (save == TRUE) {
      boot$boot <- list(predicted = predicted)
    } else {
      boot$boot <- NULL
    }
    ty <- apply(predicted$y, 2, quantile, 
                probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), 
                na.rm = TRUE)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx, y = ty)
  } else { # The following else is if resamp=0
    boot <- NULL
    boot.summary <- NULL
  }
  res <- list(real = real, boot = boot, max.distance = maxdist, 
              call = deparse(match.call()))
  class(res) <- "spline.correlog"
  res
}

#' @title Plots spline correlograms
#' @description `plot' method for class "spline.correlog".
#' @param x an object of class "spline.correlog", usually, as a result of a call to \code{spline.correlog}.
#' @param xmax the maximal distance to be plotted on the x-axis. If set to zero the maximum distance in the data will be used.
#' @param ylim limits for the y-axis (default: -1, 1).
#' @param \dots other arguments
#' @return A plot of the spline correlogram function against distance is produced. 95\% pointwise confidence (or null) envelopes are superimposed (if available).
#' @seealso \code{\link{spline.correlog}}, \code{\link{summary.spline.correlog}}
#' @keywords smooth regression
#' @export
################################################################################
plot.spline.correlog <- function(x, xmax = 0, ylim = c(-1, 1), ...) {
  ##############################################################################
  # this is the generic plot function for spline.correlog objects
  ##############################################################################
  args.default <- list(xlab = "Distance", ylab = "Correlation")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  xmax <- ifelse(xmax == 0, x$max.distance, xmax)
  do.call(plot, c(list(x = x$real$predicted$x, y = x$real$predicted$y, 
                       xlim = c(0, xmax), ylim = ylim, type = "l"), args))
  if (!is.null(x$boot$boot.summary)) {
    polygon(c(x$boot$boot.summary$predicted$x, rev(x$boot$boot.summary$predicted$x)), 
            c(x$boot$boot.summary$predicted$y["0.025", ], 
              rev(x$boot$boot.summary$predicted$y["0.975", ])), 
            col = gray(0.8), lty = 0)
  }
  lines(x$real$predicted$x, x$real$predicted$y)
  lines(c(0, max(x$real$predicted$x)), c(0, 0))
}

#' @title Summarizing spline correlograms
#' @description `summary' method for class "spline.correlog".
#' @param object an object of class "spline.correlog", usually, as a result of a call to \code{\link{spline.correlog}}.
#' @param \dots other arguments
#' @return A list summarizing spline correlograms is returned. 
#' \item{estimates}{a vector of benchmark statistics:}
#' \item{$x}{is the lowest value at which the function is = 0. If correlation is initially negative, the distance calculated appears as a negative measure.}
#' \item{$e}{is the lowest value at which the function is <= 1/e.}
#' \item{$y}{is the extrapolated value at x=0.} 
#' \item{quantiles}{a matrix summarizing the quantiles in the bootstrap (or null) distributions of the benchmark statistics.}
#' @seealso \code{\link{spline.correlog}}, \code{\link{plot.spline.correlog}}
#' @keywords smooth regression
#' @export
################################################################################
summary.spline.correlog <- function(object, ...) {
  ##############################################################################
  # this is the generic summary function for spline.correlog objects
  ##############################################################################
  xy <- cbind(object$real$x.intercept, object$real$e.intercept, object$real$y.intercept)
  dimnames(xy) <- list(c("estimate"), c("x", "e", "y"))
  yd <- apply(object$boot$boot.summary$y.intercept, 2, quantile, 
              probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
  xd <- apply(object$boot$boot.summary$x.intercept, 2, quantile, 
              probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
  ed <- apply(object$boot$boot.summary$e.intercept, 2, quantile, 
              probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
  xyd <- cbind(xd, ed, yd)
  dimnames(xyd) <- list(c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), c("x", "e", "y"))
  res <- list(call = object$call, estimate = xy, quantiles = xyd)
  res
}

#' @title Print function for spline.correlog objects
#' @description `print' method for class "spline.correlog".
#' @param x an object of class "spline.correlog", usually, as a result of a call to \code{spline.correlog} or related).
#' @param \dots other arguments
#' @return The function-call is printed to screen.
#' @seealso \code{\link{spline.correlog}}
#' @export
################################################################################
print.spline.correlog <- function(x, ...) {
  ##############################################################################
  cat("This is an object of class spline.correlog produced by the call:\n\n", x$call, 
      "\n\n Use summary() or plot() for inspection,  (or print.default() to see all the gory details).")
}
