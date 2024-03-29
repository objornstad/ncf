% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline-correlog2D.R
\name{spline.correlog2D}
\alias{spline.correlog2D}
\title{Anisotropic nonparametric (cross-)correlation function for univariate spatial data}
\usage{
spline.correlog2D(
  x,
  y,
  z,
  w = NULL,
  df = NULL,
  type = "boot",
  resamp = 1000,
  npoints = 300,
  save = FALSE,
  max.it = 25,
  xmax = FALSE,
  na.rm = FALSE,
  jitter = FALSE,
  quiet = FALSE,
  angle = c(0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5)
)
}
\arguments{
\item{x}{vector of length n representing the x coordinates.}

\item{y}{vector of length n representing the y coordinates.}

\item{z}{vector of length n representing the observation at each location.}

\item{w}{an optional second vector of length n for variable 2 (to estimate spatial or lagged cross-correlation functions).}

\item{df}{degrees-of-freedom for the spline. Default is sqrt(n).}

\item{type}{takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator}

\item{resamp}{the number of resamples for the bootstrap or the null distribution.}

\item{npoints}{the number of points at which to save the value for the spline function (and confidence envelope / null distribution).}

\item{save}{If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).}

\item{max.it}{the maximum iteration for the Newton method used to estimate the intercepts.}

\item{xmax}{If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.}

\item{na.rm}{If TRUE, NA's will be dealt with through pairwise deletion of missing values for each pair of time series -- it will dump if any one pair has less than two (temporally) overlapping observations.}

\item{jitter}{If TRUE, jitters the distance matrix to avoid problems associated with fitting the function to data on regular grids.}

\item{quiet}{If TRUE, the counter is suppressed during execution.}

\item{angle}{specifies number of cardinal directions and angles for which to calculate correlation functions. Default are 8 directions between 0 and 180.}
}
\value{
An object of class "Sncf2D" is returned. See \code{\link{Sncf2D}} for details.
}
\description{
\code{spline.correlog2D} is the function to estimate the anisotropic nonparametric correlation function in 8 (or arbitrary) directions (North - Southeast) for univariate data. Correlation functions are calculated for each different bearing. The function assumes univariate observations at each location. (use \code{\link{Sncf2D}} otherwise).
}
\details{
see \code{\link{Sncf2D}}
}
\note{
The function to estimate the UNIvariate anisotropic nonparametric (cross-)correlation function in arbitrary directions. In particular it was developed to calculate the univariate lagged cross-correlation function used in (Humston et al. 2005). Note that this 2D spline correlogram does the anisotropic analysis NOT by doing the angle-with-tolerance-wedge-style of Oden and Sokal (1986) but by projecting the the spatial coordinates of all locations on a sequence of cardinal angles (a la Sncf2D). Hence, all data points are used every time, it is only their relative distances that are changed. For example \{0, 0\} and \{0, 10\} are distance zero in the zero-degree direction but at distance 10 in the 90-degree direction.
}
\references{
Oden, N.L. and Sokal, R.R. 1986. Directional autocorrelation: an extension of spatial correlograms to two dimensions. Systematic Zoology 35: 608-617. <doi:10.2307/2413120>
  @references Humston, R., Mortensen, D. and Bjornstad, O.N. 2005. Anthropogenic forcing on the spatial dynamics of an agricultural weed: the case of the common sunflower. Journal of Applied Ecology 42: 863-872. <doi:10.1111/j.1365-2664.2005.01066.x>
}
\seealso{
\code{\link{Sncf2D}}
}
\keyword{regression}
\keyword{smooth}
