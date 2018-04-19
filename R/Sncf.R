#' @title Nonparametric (cross-)correlation function for spatio-temporal data
#' @description \code{Sncf} is the function to estimate the nonparametric (cross-)correlation function using a smoothing spline as an equivalent kernel. The function requires multiple observations at each location (use \code{\link{spline.correlog}} otherwise).
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z matrix of dimension n x p representing p observation at each location.
#' @param w an optional second matrix of dimension n x p for species 2 (to estimate the spatial cross-correlation function).
#' @param df degrees of freedom for the spline. Default is sqrt(n).
#' @param type takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param npoints the number of points at which to save the value for the spline function (and confidence envelope / null distribution).
#' @param save If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).
#' @param filter If TRUE, the Fourier filter method of Hall and coworkers is applied to ensure positive semidefiniteness of the estimator. (more work may be needed on this.)
#' @param fw If filter is TRUE, it may be useful to truncate the function at some distance w sets the truncation distance. when set to zero no truncation is done.
#' @param max.it the maximum iteration for the Newton method used to estimate the intercepts.
#' @param xmax If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values for each pair of time series -- it will dump if any one pair has less than two (temporally) overlapping observations.
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param circ If TRUE, the observations are assumed to be angular (in radians), and circular correlation is used.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "Sncf" is returned, consisting of the following components: 
#' \item{real}{the list of estimates from the data.}
#' \item{$cbar}{the regional average correlation.}
#' \item{$x.intercept}{the lowest value at which the function is = 0. If correlation is initially negative, the distance is given as negative.}
#' \item{$e.intercept}{the lowest value at which the function 1/e.}
#' \item{$y.intercept}{the extrapolated value at x=0 (nugget).}
#' \item{$cbar.intercept}{distance at which regional average correlation is reach.}
#' \item{$predicted$x}{the x-axes for the fitted covariance function.}
#' \item{$predcited$y}{the values for the covariance function.}
#' \item{boot}{a list with the analogous output from the bootstrap or null distribution.} 
#' \item{$summary}{gives the full vector of output for the x.intercept, y.intercept, e.intercept, cbar.intercept, cbar and a quantile summary for the resampling distribution.}
#' \item{$boot}{If save=TRUE, the full raw matrices from the resampling is saved.}
#' \item{max.distance}{the maximum spatial distance considered.}
#' @details Missing values are allowed -- values are assumed missing at random. 
#' 
#'   The circ argument computes a circular version of the Pearson's product moment correlation (see \code{\link{cor2}}). This option is to calculate the 'nonparametric phase coherence function' (Grenfell et al. 2001)
#' @references Hall, P. & Patil, P. (1994) Properties of nonparametric estimators of autocovariance for stationary random fields. Probability Theory and Related Fields, 99:399-424. \url{https://doi.org/10.1007/BF01199899}
#' 
#'   Hall, P., Fisher, N.I. & Hoffmann, B. (1994) On the nonparametric estimation of covariance functions. Annals of Statistics, 22:2115-2134. \url{https://doi.org/10.1214/aos/1176325774}
#'   
#'   Bjornstad, O.N. & Falck, W. (2001) Nonparametric spatial covariance functions: estimation and testing. Environmental and Ecological Statistics, 8:53-70. \url{https://doi.org/10.1023/A:1009601932481}
#'   
#'   Bjornstad, O.N., Ims, R.A. & Lambin, X. (1999) Spatial population dynamics: Analysing patterns and processes of population synchrony. Trends in Ecology and Evolution, 11:427-431. \url{https://doi.org/10.1016/S0169-5347(99)01677-8}
#'   
#'   Bjornstad, O. N., and J. Bascompte. (2001) Synchrony and second order spatial correlation in host-parasitoid systems. Journal of Animal Ecology 70:924-933. \url{https://doi.org/10.1046/j.0021-8790.2001.00560.x}
#'   
#'   Grenfell, B.T., Bjornstad, O.N., & Kappey, J. (2001) Travelling waves and spatial hierarchies in measles epidemics. Nature 414:716-723. \url{https://doi.org/10.1038/414716a}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{summary.Sncf}}, \code{\link{plot.Sncf}}, \code{\link{Sncf2D}}, \code{\link{Sncf.srf}}
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[, 2]
#' # z data from an exponential random field
#' z <- cbind(
#'    rmvn.spa(x = x, y = y, p = 2, method = "exp"),
#'    rmvn.spa(x = x, y = y, p = 2, method = "exp")
#'   )
#' # w data from a gaussian random field
#' w <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "gaus"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "gaus")
#'   )
#' # multivariate nonparametric covariance function
#' fit1 <- Sncf(x = x, y = y, z = z, resamp = 0)
#' \dontrun{plot.Sncf(fit1)}
#' summary(fit1)
#' 
#' # multivariate nonparametric cross-covariance function
#' fit2 <- Sncf(x = x, y = y, z = z, w = w, resamp = 0)
#' \dontrun{plot(fit2)}
#' summary(fit2)
#' @keywords smooth regression spatial
#' @export
################################################################################
Sncf <- function(x, y, z, w = NULL, df = NULL, type = "boot", resamp = 1000, 
                 npoints = 300, save = FALSE, filter = FALSE, fw = 0, max.it = 25, 
                 xmax = FALSE, na.rm = FALSE, latlon = FALSE, circ = FALSE, 
                 quiet = FALSE) {
  ##############################################################################
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
  
  if (is.null(w)) {
    # This generates correlations
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z) + 0
    moran <- cor2(t(z), circ = circ)
  } else {
    # This generates cross-correlations
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z) + 0
    w <- as.matrix(w) + 0
    moran <- cor2(t(z), t(w), circ = circ)
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
  sel <- is.finite(v) & is.finite(u)
  u <- u[sel]
  v <- v[sel]
  v <- v[u <= maxdist]
  u <- u[u <= maxdist]
  
  xpoints  <- seq(0, maxdist, length = npoints)
  out <- gather(u = u, v = v, w = w, moran = moran, df = df, xpoints = xpoints, 
                filter = filter, fw = fw)
  # End of spline fit
  real <- list(cbar = out$cbar, x.intercept = out$xint, e.intercept = out$eint, 
               y.intercept = out$yint, cbar.intercept = out$cint,
               predicted = list(x = matrix(out$x, nrow = 1),
                                y = matrix(out$y, nrow = 1)))
  
  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if (resamp != 0) {
    # here is the bootstrapping/randomization
    boot$boot.summary$x.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$cbar.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$cbar <- matrix(NA, nrow = resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), 
                      y = matrix(NA, nrow = resamp, ncol = npoints))
    type <- charmatch(type, c("boot", "perm"), nomatch = NA)
    if (is.na(type))
      stop("method should be \"boot\", or \"perm\"")
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (quiet & any(i == whn))	{
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      if (type == 1) {
        trekkx <- sample(1:n, replace = TRUE)
        trekky <- trekkx
      }
      
      if (type == 2) {
        trekky <- sample(1:n, replace = FALSE)
        trekkx <- 1:n
      }
      
      xdistb <- xdist[trekkx, trekkx]
      
      if (is.null(w)) {
        triang <- lower.tri(xdistb)
      } else {
        triang <- is.finite(xdistb)
      }
      
      xdistb <- xdistb[triang]
      moranb <- moran[trekky, trekky][triang]
      
      if (type == 1 & is.null(w)) {
        moranb <- moranb[!(xdistb == 0)]
        xdistb <- xdistb[!(xdistb == 0)]
      }
      
      u <- xdistb
      v <- moranb
      sel <- is.finite(v) & is.finite(u)
      u <- u[sel]
      v <- v[sel]
      v <- v[u <= maxdist]
      u <- u[u <= maxdist]
      
      out <- gather(u = u, v = v, w = w, moran = moranb, df = df, xpoints = xpoints, 
                    filter = filter, fw = fw)
      
      boot$boot.summary$cbar[i, 1] <- out$cbar
      boot$boot.summary$y.intercept[i, 1] <- out$yint
      boot$boot.summary$x.intercept[i, 1] <- out$xint
      boot$boot.summary$e.intercept[i, 1] <- out$eint
      boot$boot.summary$cbar.intercept[i, 1] <- out$cint
      predicted$x[1, ] <- out$x
      predicted$y[i, ] <- out$y
      
    }
    # end of bootstrap loop!
    
    if (save == TRUE) {
      boot$boot <- list(predicted = predicted)
    } else {
      boot$boot <- NULL
    }
    
    ty <- apply(predicted$y, 2, quantile, probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5,
                                                    0.75, 0.9, 0.95, 0.975, 1), 
                na.rm = TRUE)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), 
                         NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx,y = ty)
  } else { # The following else is if resamp = 0
    boot <- NULL
    boot.summary <- NULL
  }
  
  res <- list(real = real, boot = boot, max.distance = maxdist, 
              call = deparse(match.call()))
  class(res) <- "Sncf"
  res
}

#' @title Plots nonparametric spatial correlation-functions
#' @description `plot' method for class "Sncf".
#' @param x an object of class "Sncf", usually, as a result of a call to \code{Sncf} (or \code{Sncf.srf}).
#' @param xmax the maximal distance to be plotted on the x-axis. If set to zero the maximum distance in the data will be used.
#' @param ylim limits for the y-axis (default: -1, 1).
#' @param add If TRUE the plot is added on to the previous graph.
#' @param \dots other arguments
#' @return A plot of the nonparametric spatial covariance function (with CI's if boostrapps are available)
#' @seealso \code{\link{Sncf}}, \code{\link{plot.Sncf}}, \code{\link{Sncf.srf}}, \code{\link{summary.Sncf}}
#' @keywords smooth regression
#' @export
################################################################################
plot.Sncf <- function(x, xmax = 0, ylim = c(-1, 1), add = FALSE, ...) {
  ##############################################################################
  args.default <- list(xlab = "Distance", ylab = "Correlation")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  xmax <- ifelse(xmax == 0, x$max.distance, xmax)
  cbar <- x$real$cbar
  if (!add) {
    do.call(plot, c(list(x = x$real$predicted$x, y = x$real$predicted$y, 
                         xlim = c(0, xmax), ylim = ylim, type = "l"), args))
  }
  if (!is.null(x$boot$boot.summary)) {
    polygon(c(x$boot$boot.summary$predicted$x, rev(x$boot$boot.summary$predicted$x)), 
            c(x$boot$boot.summary$predicted$y["0.025", ], 
              rev(x$boot$boot.summary$predicted$y["0.975", ])), col = gray(0.8), 
            lty = 0)
  }
  lines(x$real$predicted$x, x$real$predicted$y)
  lines(c(0, max(x$real$predicted$x)), c(0, 0))
  lines(c(0, max(x$real$predicted$x)), c(cbar, cbar))
}

#' @title Print function for Sncf objects
#' @description `print' method for class "Sncf".
#' @param x an object of class "Sncf", usually, as a result of a call to \code{Sncf} or related).
#' @param \dots other arguments
#' @return The function-call is printed to screen.
#' @seealso \code{\link{Sncf}}
#' @export
################################################################################
print.Sncf <- function(x, ...) {
  ##############################################################################
  cat("This is an object of class Sncf produced by the call:\n\n", x$call, 
      "\n\n Use summary() or plot() for inspection,  
      (or print.default() to see all the gory details).")
}

#' @title Summarizing nonparametric spatial correlation-functions
#' @description `summary' method for class "Sncf".
#' @param object an object of class "Sncf", usually, as a result of a call to \code{\link{Sncf}} (or \code{\link{Sncf.srf}}).
#' @param \dots other arguments
#' @return A list summarizing the nonparametric (cross-)covariance function is returned. 
#' \item{Regional.synch}{the regional mean (cross-)correlation.}
#' \item{Squantile}{the quantile distribution from the resampling for the regional correlation.}
#' \item{estimates}{a vector of benchmark statistics:}
#' \item{$x}{is the lowest value at which the function is = 0. If correlation is initially negative, the distance calculated appears as a negative measure.}
#' \item{$e}{is the lowest value at which the function is <= 1/e.}
#' \item{$y}{is the extrapolated value at x=0.} 
#' \item{$cbar}{is the shortest distance at which function is = regional mean correlation.}
#' \item{quantiles}{a matrix summarizing the quantiles in the bootstrap (or null) distributions of the benchmark statistics.}
#' @seealso \code{\link{Sncf}}, \code{\link{plot.Sncf}}
#' @keywords smooth regression
#' @export
################################################################################
summary.Sncf <- function(object, ...) {
  ##############################################################################
  # this is the generic summary function for Sncf objects
  ##############################################################################
  xy <- cbind(object$real$x.intercept, object$real$e.intercept, 
              object$real$y.intercept, object$real$cbar.intercept)
  dimnames(xy) <- list(c("intercepts"), c("x", "e","y", "cbar"))
  if (!is.null(object$boot$boot.summary)) {
    yd <- apply(object$boot$boot.summary$y.intercept, 2, quantile, 
                probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
    xd <- apply(object$boot$boot.summary$x.intercept, 2, quantile, 
                probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
    ed <- apply(object$boot$boot.summary$e.intercept, 2, quantile, 
                probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
    synchd <- quantile(object$boot$boot.summary$cbar[, 1], 
                       probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
    cbard <- quantile(object$boot$boot.summary$cbar.intercept[, 1], 
                      probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE)
    xyd <- cbind(xd, ed, yd, cbard)
    dimnames(xyd) <- list(c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), c("x", "e", "y", "cbar"))
  }
  
  if (is.null(object$boot$boot.summary)) {
    synchd <- NULL
    xyd <- NULL	
  }
  
  res <- list(call = object$call, Regional.synch = object$real$cbar, Squantile = synchd, 
              estimates = xy, quantiles = xyd)
  res
}

#' @title Nonparametric (Cross-)Covariance Function from stationary random fields
#' @description \code{Sncf.srf} is the function to estimate the nonparametric for spatio-temporal data from fully stationary random fields (i.e. marginal expectation and variance identical for all locations; use \code{\link{Sncf}} otherwise).
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z matrix of dimension n x p representing p observation at each location.
#' @param w an optional second matrix of dimension n x p for variable 2 (to estimate the spatial cross-correlation function).
#' @param avg supplies the marginal expectation of the Markov random field; if TRUE, the sample mean (across the markovian field) is used.
#' @param avg2 optionally supplies the marginal expectation of the Markov random field for optional variable 2; if TRUE, the sample mean is used.
#' @param corr If TRUE, the covariance function is standardized by the marginal variance (across the markovian field) to return a correlation function (alternatively the covariance function is returned).
#' @param df degrees of freedom for the spline. Default is sqrt(n).
#' @param type takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param npoints the number of points at which to save the value for the spline function (and confidence envelope / null distribution).
#' @param save If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).
#' @param filter If TRUE, the Fourier filter method of Hall and coworkers is applied to ensure positive semidefiniteness of the estimator. (more work may be needed on this.)
#' @param fw If filter is TRUE, it may be useful to truncate the function at some distance w sets the truncation distance. When set to zero no truncation is done.
#' @param max.it the maximum iteration for the Newton method used to estimate the intercepts.
#' @param xmax If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.
#' @param jitter If TRUE, jitters the distance matrix, to avoid problems associated with fitting the function to data on regular grids.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "Sncf" (or "Sncf.cov") is returned. See \code{\link{Sncf}} for details.
#' @details If \code{corr = F}, an object of class "Sncf.cov" is returned. Otherwise the class is "Sncf".
#' 
#'   \code{Sncf.srf} is a function to estimate the nonparametric (cross-)covariance function (as discussed in Bjornstad and Bascompte 2001) for data from a fully stationary random fields. I have found it useful to estimate the (cross-)covariance functions in synthetic data.
#' @references Bjornstad, O. N., and J. Bascompte. (2001) Synchrony and second order spatial correlation in host-parasitoid systems. Journal of Animal Ecology 70:924-933. \url{https://doi.org/10.1046/j.0021-8790.2001.00560.x}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{Sncf}}, \code{\link{summary.Sncf}}, \code{\link{plot.Sncf}}, \code{\link{plot.Sncf.cov}}
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
#' # multivariate nonparametric covariance function
#' fit1 <- Sncf.srf(x = x, y = y, z = z, avg = NULL, corr = TRUE, resamp = 0) 
#' \dontrun{plot(fit1)} 
#' summary(fit1)
#' 
#' # multivariate nonparametric cross-covariance function (with known
#' # marginal expectation of zero for both z and w
#' fit2 <- Sncf.srf(x = x, y = y, z = z, w = w, avg = 0, avg2 = 0, corr = FALSE, 
#'                  resamp = 0)
#' \dontrun{plot(fit2)} 
#' summary(fit2)
#' @keywords smooth regression
#' @export
################################################################################
Sncf.srf <- function(x, y, z, w = NULL, avg = NULL, avg2 = NULL, corr = TRUE, 
                     df = NULL, type = "boot", resamp = 0, npoints = 300, 
                     save = FALSE, filter = FALSE, fw = 0, max.it = 25, 
                     xmax = FALSE, jitter = FALSE, quiet = FALSE) {
  ##############################################################################
  # Sncf.srf is the function to estimate the nonparametric covariance function for a 
  # stationary random field (expectation and variance identical). The function uses a
  # smoothing spline as an equivalent kernel) as discussed in 
  # Bjornstad et al. (1999; Trends in Ecology and Evolution 14:427-431)
  ##############################################################################
  p <- dim(z)[2]
  n <- dim(z)[1]
  
  if (is.null(df)) {
    df <- sqrt(n)
  }
  
  if (is.null(avg)) {
    avg <- mean(as.vector(z), na.rm = TRUE)
    
    if (!is.null(w)) {
      avg2 <- mean(as.vector(w), na.rm = TRUE)
    }
  }
  
  sca <- 1
  sca2 <- 1
  if (corr == TRUE) {
    sca <- sqrt(var(as.vector(z)))
    
    if (!is.null(w)) {
      sca2 <- sqrt(var(as.vector(w)))
    }
  }
  
  # generates distance matrices
  xdist <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
  
  if (jitter == TRUE) {
    # xdist <- jitter(xdist)
    xdist <- apply(xdist, 2, jitter)
  }
  
  if (is.null(w)) {
    moran <- crossprod((t(z) - avg)/(sca))/p
  } else {
    moran <- crossprod((t(z) - avg)/(sca), (t(w) - avg2)/(sca2))/p
  }
  
  maxdist <- ifelse(!xmax, max(xdist), xmax)
  
  # Spline fit
  
  # The spline function
  if (is.null(w)) {
    triang <- lower.tri(xdist)
  } else {
    triang <- is.finite(xdist)
  }
  
  u <- xdist[triang]
  v <- moran[triang]
  v <- v[u <= maxdist]	
  u <- u[u <= maxdist]
  
  xpoints <- seq(0, maxdist, length = npoints)
  out <- gather(u = u, v = v, w = w, moran = moran, df = df, xpoints = xpoints, 
                filter = filter, fw = fw)
  # End of spline fit
  real <- list(cbar = out$cbar, x.intercept = out$xint, e.intercept = out$eint, 
               y.intercept = out$yint, cbar.intercept = out$cint,
               predicted = list(x = matrix(out$x, nrow = 1),
                                y = matrix(out$y, nrow = 1)))
  
  boot <- list(NULL)
  boot$boot.summary <- list(NULL)
  if (resamp != 0) {
    # here is the bootstrapping/randomization
    boot$boot.summary$x.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$cbar.intercept <- matrix(NA, nrow = resamp, ncol = 1)
    boot$boot.summary$cbar <- matrix(NA, nrow = resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), 
                      y = matrix(NA, nrow = resamp, ncol = npoints))
    type <- charmatch(type, c("boot", "perm"), nomatch = NA)
    if (is.na(type))
      stop("method should be \"boot\", or \"perm\"")
    for (i in 1:resamp) {
      whn <- pretty(c(1, resamp), n = 10)
      if (quiet & any(i == whn))	{
        cat(i, " of ", resamp, "\r")
        flush.console()
      }
      
      if (type == 1) {
        trekkx <- sample(1:n, replace = TRUE)
        trekky <- trekkx
      }
      
      if (type == 2) {
        trekky <- sample(1:n, replace = FALSE)
        trekkx <- 1:n
      }
      
      xdistb <- xdist[trekkx, trekkx]
      
      if (is.null(w)) {
        triang <- lower.tri(xdistb)
      } else {
        triang <- is.finite(xdistb)
      }
      
      xdistb <- xdistb[triang]
      moranb <- moran[trekky, trekky][triang]
      
      if (type == 1 & is.null(w)) {
        moranb <- moranb[!(xdistb == 0)]
        xdistb <- xdistb[!(xdistb == 0)]
      }
      
      u <- xdistb
      v <- moranb
      v <- v[u <= maxdist]	
      u <- u[u <= maxdist]
      out <- gather(u = u, v = v, w = w, moran = moranb, df = df, xpoints = xpoints, 
                    filter = filter, fw = fw)
      
      boot$boot.summary$cbar[i, 1]  <- out$cbar
      boot$boot.summary$y.intercept[i, 1] <- out$yint
      boot$boot.summary$x.intercept[i, 1] <- out$xint
      boot$boot.summary$e.intercept[i, 1] <- out$eint
      boot$boot.summary$cbar.intercept[i, 1] <- out$cint
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
  if (corr) {
    class(res) <- "Sncf"
  } else {
    class(res) <- "Sncf.cov"
  }
  
  res
}

#' @title Plots nonparametric spatial covariance-functions
#' @description `plot' method for class "Sncf.cov".
#' @param x an object of class "Sncf.cov", usually, as a result of a call to \code{Sncf.srf} (with \code{corr} = FALSE).
#' @param xmax the maximal distance to be plotted on the x-axis. If set to zero the maximum distance in the data will be used.
#' @param \dots other arguments
#' @return A plot of the nonparametric spatial covariance function (with CI's if boostrapps are available)
#' @seealso \code{\link{Sncf.srf}}, \code{\link{plot.Sncf}}
#' @keywords smooth regression
#' @export
################################################################################
plot.Sncf.cov <- function(x, xmax = 0, ...) {
  ##############################################################################
  args.default <- list(xlab = "Distance", ylab = "Covariance")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  xmax <- ifelse(xmax == 0, max(x$real$predicted$x), xmax)
  do.call(plot, c(list(x = x$real$predicted$x, y = x$real$predicted$y, 
                       xlim = c(0, xmax), type = "l"), args))
  if (!is.null(x$boot$boot.summary)) {
    polygon(c(x$boot$boot.summary$predicted$x, rev(x$boot$boot.summary$predicted$x)), 
            c(x$boot$boot.summary$predicted$y["0.025", ], 
              rev(x$boot$boot.summary$predicted$y["0.975", ])), col = gray(0.8), 
            lty = 0)
  }
  lines(x$real$predicted$x, x$real$predicted$y)
  lines(c(0, max(x$real$predicted$x)), c(0, 0))
}
