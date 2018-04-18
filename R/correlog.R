#' @title Uni- and multivariate spatial correlograms
#' @description \code{correlog} is the function to estimate spatial (cross-)correlograms. Either univariate or multivariate (time seres) for each site can be used.
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z vector of length n or matrix of dimension n x p representing p observation at each location.
#' @param w an optional second variable with idenitical dimension to z (to estimate cross-correlograms).
#' @param increment increment for the uniformly distributed distance classes.
#' @param resamp the number of permutations under the null to assess level of significance.
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "correlog" is returned, consisting of the following components: 
#' \item{correlation}{the value for the moran (or Mantel) similarity.}
#' \item{mean.of.class}{the actual average of the distances within each distance class.}
#' \item{nlok}{the number of pairs within each distance class.}
#' \item{x.intercept}{the interpolate x.intercept of Epperson (1993).}
#' \item{p}{the permutation two-sided p-value for each distance-class.}
#' \item{corr0}{If a cross-correlogram is calculated, corr0 gives the empirical cross-correlation at distance zero.}
#' @details The spatial (cross-)correlogram and Mantel (cross-)correlogram estimates the spatial dependence at discrete distance classes. 
#' 
#'  The regionwide similarity forms the reference line (the zero-line); the x-intercept is thus the distance at which object are no more similar than that expected by-chance-alone across the region.
#'  
#'  If the data are univariate, the spatial dependence is measured by Moran's \emph{I}. If it is multivariate, it is measured by the \emph{centred} Mantel statistic. (Use \code{\link{correlog.nc}} if the non-centered multivariate correlogram is desired).
#'  
#'  Missing values are allowed -- values are assumed missing at random.
#' @references Bjornstad, O.N., Ims, R.A. & Lambin, X. (1999) Spatial population dynamics: Analysing patterns and processes of population synchrony. Trends in Ecology and Evolution, 11, 427-431. \url{https://doi.org/10.1016/S0169-5347(99)01677-8}
#' 
#'   Bjornstad, O.N. & Falck, W. (2001) Nonparametric spatial covariance functions: estimation and testing. Environmental and Ecological Statistics, 8:53-70. \url{https://doi.org/10.1023/A:1009601932481}
#'   
#'   Epperson, B.K. (1993) Recent advances in correlation studies of spatial patterns of genetic variation. Evolutionary Biology, 27, 95-155. \url{https/doi.org/10.1007/978-1-4615-2878-4_4}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{plot.correlog}}, \code{\link{spline.correlog}}, \code{\link{correlog.nc}}
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
#' # Spatial correlogram 
#' fit1 <- correlog(x = x, y = y, z = z[, 1], increment = 2, resamp = 0) 
#' \dontrun{plot(fit1)}
#' 
#' # Mantel correlogram 
#' fit2 <- correlog(x = x, y = y, z = z, increment = 2, resamp = 0) 
#' \dontrun{plot(fit2)}
#' 
#' # Mantel cross-correlogram 
#' fit3 <- correlog(x = x, y = y, z = z, w = w, increment = 2, resamp = 0) 
#' \dontrun{plot(fit3)}
#' @keywords spatial
#' @export
#' @importFrom grDevices gray
#' @importFrom graphics lines par plot points polygon symbols text title
#' @importFrom stats cor fft lm na.omit predict quantile rnorm sd smooth.spline uniroot var
#' @importFrom utils flush.console
##############################################################################################
correlog<-function(x, y, z, w=NULL, increment, resamp = 1000, latlon = FALSE, na.rm = FALSE, quiet=FALSE){
  ##############################################################################################
  #correlog estimates the spatial correlogram (if z is univariate)
  #or the Mantel correlogram (if z is multivariate), or the (uni-/multivariate)
  #cross-correlogram (if the optional w data set is given).
  #######################################################################################
  
  NAO <- FALSE
  
  #check for missing values
  if(any(!is.finite(unlist(z)))) {
    if(na.rm){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  
  
  multivar <- !is.null(dim(z))		#test whether z is univariate or multivariate
  
  if(multivar == TRUE) {
    warning("Response is multivariate: the correlation matrix will be centered on zero. Use correlog.nc() for the non-centered correlogram")
    n <- length(z[, 1])
    p <- length(z[1,  ])
    z <- as.matrix(z)+0
    
    if(is.null(w)){
      moran <- cor2(t(z), circ=FALSE)
      moran <- moran[lower.tri(moran)]
      moran <- moran - mean(moran, na.rm= TRUE)
    }
    
    else {
      w <- as.matrix(w)+0
      moran <- cor2(t(z), t(w), circ=FALSE)
      zero <- mean(diag(moran), na.rm= TRUE)
      moran <- moran[row(moran)!=col(moran)]
      moran <- moran - mean(moran, na.rm= TRUE)
    }
  }
  
  else {
    n <- length(z)
    z <- as.vector(z)+0
    zscal <- (scale(z, center = TRUE, scale = TRUE)[, 1])/(sqrt((n-1)/n))
    
    if(is.null(w)){
      moran <- t(outer(zscal, zscal))
      moran <- moran[lower.tri(moran)]
    }
    
    else {
      wscal <- (scale(w, center = TRUE, scale = TRUE)[, 1])/(sqrt((n-1)/n))
      zw <- c(zscal,wscal)
      #this may be a slight hack
      moran <- t(outer(zw, zw))[1:n,(n+1):(2*n)]
      zero <- mean(diag(moran), na.rm= TRUE)
      moran <- moran[row(moran)!=col(moran)]
    }
  }
  
  
  #then generating geographic distances
  if(latlon){
    #these are geographic distances from lat-lon coordinates
    dmat <- gcdist(x,y)
  }
  
  else{
    #these are geographic distances from euclidian coordinates
    dmat <- sqrt(outer(x,x, "-")^2+outer(y,y,"-")^2)
  }
  
  if(resamp != 0){
    dmat2 <- dmat
    moran2 <- moran
  }
  
  if(is.null(w)){
    dmat <- dmat[lower.tri(dmat)]
  }
  else {
    dmat <- dmat[row(dmat)!=col(dmat)]
  }
  
  dkl <- ceiling(dmat/increment) 	#generates the distance matrices
  nlok <- sapply(split(moran, dkl), length)
  dmean <- sapply(split(dmat, dkl), mean, na.rm = TRUE)
  moran <- sapply(split(moran, dkl), mean, na.rm = TRUE)
  ly <- 1:length(dmean)
  x <- c(dmean[ly[moran < 0][1]], dmean[ly[moran < 0][1] - 1])
  y <- c(moran[ly[moran < 0][1] - 1], moran[ly[moran < 0][1]])
  if(moran[1] < 0) {
    tmp <- 0
  }
  else {
    tmp <- lm(x ~ y)[[1]][1]
  }
  
  p<-NULL
  
  if(resamp != 0){
    
    perm <- matrix(NA, ncol = length(moran), nrow = resamp)
    
    for(i in 1:resamp){
      whn=pretty(c(1,resamp), n=10)
      if(! quiet & any(i==whn)){
        cat(i, " of ", resamp, "\r")
        flush.console()}
      
      
      trekk <-sample(1:n)
      dma <- dmat2[trekk,trekk]
      mor <- moran2
      
      if(is.null(w)){
        dma <- dma[lower.tri(dma)]
      }
      else {
        dma <- dma[row(dma)!=col(dma)]
      }
      
      dkl <- ceiling(dma/increment)	#generates the distance matrices
      perm[i,] <- sapply(split(mor, dkl), mean, na.rm = TRUE)
    }
    
    p=(apply(moran<=t(perm),1,sum))/(resamp+1)
    p=apply(cbind(p, 1-p), 1, min) + 1/(resamp+1)
  }
  
  res <- list(n = nlok, mean.of.class = dmean, correlation = moran,
              x.intercept = tmp, p=p, call=deparse(match.call()))
  
  if(!is.null(w)){
    res$corr0 <- zero
  }
  class(res) <- "correlog"
  res
}

#' @title Plots spatial correlograms
#' @description `plot' method for class "correlog".
#' @param x an object of class "correlog", usually, as a result of a call to \code{correlog} or \code{correlog.nc}.
#' @param \dots other arguments
#' @return A spatial or Mantel (cross-correlogram) is plotted. 
#' 
#'   If a permutation test was performed, values significant at a nominal (two-sided) 5\%-level will be respresented by filled circles and non-significant values by open circles.
#' @seealso \code{\link{correlog}}, \code{\link{correlog.nc}}
#' @keywords spatial
#' @export
#########################################################################################
plot.correlog<-function(x, ...){
  ##############################################################################################
  #this is the generic plot function for correlog objects
  #sigificant values are represented by filled cirles
  ##############################################################################################
  plot(x$mean.of.class, x$correlation, ylab='correlation', xlab='distance (mean-of-class)')
  lines(x$mean.of.class, x$correlation)
  if(!is.null(x$p)){
    points(x$mean.of.class[x$p<0.025], x$correlation[x$p<0.025], pch=21, bg="black")
  }
  title("Correlogram")
}

#' @title Non-cenetered spatial (cross-)correlogram
#' @description \code{correlog.nc} is the function to estimate the non-centred (cross-)correlogram. The noncentred correlogram provides estimates of the spatial correlation for discrete distance classes. The function requires multiple observations at each location (use \code{\link{correlog}} otherwise).
#' @param x vector of length n representing the x coordinates (or longitude; see latlon).
#' @param y vector of length n representing the y coordinates (or latitude).
#' @param z a matrix of dimension n x p representing p (>1) observation at each location.
#' @param w an optional second variable with idenitical dimension to z (to estimate cross-correlograms).
#' @param increment increment for the uniformly distributed distance classes.
#' @param resamp the number of permutations under the null to assess level of significance.
#' @param latlon If TRUE, coordinates are latitude and longitude.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "correlog" is returned, consisting of the following components: 
#' \item{correlation}{the value for the moran (or Mantel) similarity.}
#' \item{mean.of.class}{the actual average of the distances within each distance class.}
#' \item{nlok}{the number of pairs within each distance class.}
#' \item{x.intercept}{the interpolate x.intercept of Epperson (1993).}
#' \item{p}{the permutation p-value for each distance-class.}
#' \item{corr0}{If a cross-correlogram is calculated, corr0 gives the empirical within-patch cross-correlation.}
#' @details The non-centred correlogram estimates spatial dependence at discrete distance  classes. The method corresponds to the modified correlogram of Koenig & Knops(1998), but augumented to potentially estimate the cross-correlogram). The function requires multiple observations at each location. Missing values is allowed in the multivariate case (pairwise deletion will be used).
#'   
#'   Missing values are allowed -- values are assumed missing at random.
#' @references Bjornstad, O.N., Ims, R.A. & Lambin, X. (1999) Spatial population dynamics: Analysing patterns and processes of population synchrony. Trends in Ecology and Evolution, 11, 427-431. \url{https://doi.org/10.1016/S0169-5347(99)01677-8}
#'   
#'   Koenig, W.D. & Knops, J.M.H. (1998) Testing for spatial autocorrelation in ecological studies. Ecography, 21, 423-429. \url{https://doi.org/10.1111/j.1600-0587.1998.tb00407.x}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{plot.correlog}}, \code{\link{correlog}}
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
#' # noncentered (Mantel) correlogram 
#' fit1 <- correlog.nc(x = x, y = y, z = z, increment = 2, resamp = 500)
#' \dontrun{plot.correlog(fit1)}
#' @keywords spatial
#' @export
##############################################################################################
correlog.nc<-function(x, y, z, w=NULL, increment, resamp = 1000, na.rm = FALSE, latlon=FALSE, quiet=FALSE){
  ##############################################################################################
  #correlog.nc estimates the noncentred correlogram
  #and cross-correlogram. Bjornstad et al. (1999; Trends in Ecology and
  #Evolution 14:427-431)
  #The function requires mulitple observations at each location (use
  #correlog otherwise).
  #######################################################################################
  
  NAO <- FALSE
  
  #check for missing values
  if(any(!is.finite(unlist(z)))) {
    if(na.rm){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  
  #then generating geographic distances
  if(latlon){
    #these are geographic distances from lat-lon coordinates
    dmat <- gcdist(x,y)
  }
  
  else{
    #these are geographic distances from euclidian coordinates
    dmat <- sqrt(outer(x,x, "-")^2+outer(y,y,"-")^2)
  }
  
  if(is.null(w)){
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z)+0
    
    moran <- cor2(t(z), circ=FALSE)
  }
  
  else {
    #This generates the moran distances for cross-correlation
    #the odd adding of zero is just to ensure that all vectors
    #are treated as numeric
    n <- dim(z)[1]
    p <- dim(z)[2]
    z <- as.matrix(z)+0
    w <- as.matrix(w)+0
    
    moran <- cor2(t(z), t(w), circ=FALSE)
  }
  
  if(resamp != 0){
    dmat2 <- dmat
    #moran2 <- moran
  }
  
  if(is.null(w)){
    dmat <- dmat[lower.tri(dmat)]
    moran <- moran[lower.tri(moran)]
  }
  else {
    dmat <- dmat[row(dmat)!=col(dmat)]
    zero <- mean(diag(moran), na.rm= TRUE)
    moran <- moran[row(moran)!=col(moran)]
  }
  
  if(resamp != 0){
    #dmat2 <- dmat
    moran2 <- moran
  }
  
  dkl <- ceiling(dmat/increment) 	#generates the distance matrices
  nlok <- sapply(split(moran, dkl), length)
  dmean <- sapply(split(dmat, dkl), mean, na.rm = TRUE)
  moran <- sapply(split(moran, dkl), mean, na.rm = TRUE)
  ly <- 1:length(dmean)
  x <- c(dmean[ly[moran < 0][1]], dmean[ly[moran < 0][1] - 1])
  y <- c(moran[ly[moran < 0][1] - 1], moran[ly[moran < 0][1]])
  
  if(moran[1] < 0) {
    tmp <- 0
  }
  else {
    if(sum(moran<0)>0){
      tmp <- lm(x ~ y)[[1]][1]
    }
    else{
      tmp <- NA
    }
  }
  
  p<-NULL
  
  if(resamp != 0){
    perm <- matrix(NA, ncol = length(moran), nrow = resamp)
    
    for(i in 1:resamp){
      whn=pretty(c(1,resamp), n=10)
      if(! quiet & any(i==whn)){
        cat(i, " of ", resamp, "\r")
        flush.console()}
      
      
      
      trekk <-sample(1:n)
      dma <- dmat2[trekk,trekk]
      mor <- moran2
      
      if(is.null(w)){
        dma <- dma[lower.tri(dma)]
      }
      else {
        dma <- dma[row(dma)!=col(dma)]
      }
      
      dkl <- ceiling(dma/increment)	#generates the distance matrices
      perm[i,] <- sapply(split(mor, dkl), mean, na.rm = TRUE)
    }
    
    p=(apply(moran<=t(perm),1,sum))/(resamp+1)
    p=apply(cbind(p, 1-p), 1, min) + 1/(resamp+1)
  }
  
  res <- list(n = nlok, mean.of.class = dmean, correlation = moran,
              x.intercept = tmp, p = p, call=deparse(match.call()))
  if(!is.null(w)){
    res$corr0 <- zero
  }
  class(res) <- "correlog"
  res
}

#' @title Mantel (cross-)correlograms
#' @description \code{mantel.correlog} is the function to calculate a Mantel (cross-)correlogram. The function requires two (or three) matrices.
#' @param dmat a matrix representing distance.
#' @param zmat a matrix representing similarity.
#' @param wmat an optional third matrix of similarities to calculate a Mantel cross-correlograms.
#' @param increment increment for the uniformly distributed distance classes.
#' @param resamp the number of permutations under the null to assess level of significance.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "correlog" is returned, consisting of the following components: 
#' \item{correlation}{the value for the moran (or Mantel) similarity.}
#' \item{mean.of.class}{the actual average of the distances within each distance class.}
#' \item{nlok}{the number of pairs within each distance class.}
#' \item{x.intercept}{the interpolate x.intercept of Epperson (1993).}
#' \item{p}{the permutation two-sided p-value for each distance-class.}
#' \item{corr0}{If a cross-correlogram is calculated, corr0 gives the empirical cross-correlation at distance zero.}
#' @details The function calculates Mantel (cross-)correlograms at discrete distance classes from two (or three) matrixes. The first is the matrix of distances and the second is a matrix of similarities. The optional third matrix is an additional similarity matrix to be used to calculate a Mantel cross-correlogram.  Missing values are allowed -- values are assumed missing at random.
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{plot.correlog}}
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
#' w <- cbind(rmvn.spa(
#'   x = x, y = y, p = 2, method = "gaus"), 
#'   rmvn.spa(x = x, y = y, p = 2, method = "gaus")
#'   )
#' 
#' # Make distance and similarity matrices
#' zmat <- cor(t(z))
#' wmat <- cor(t(w))
#' dmat <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2)
#' 
#' # Mantel correlogram 
#' fit1 <- mantel.correlog(dmat = dmat, zmat = zmat, increment = 2, quiet = TRUE, 
#' resamp = 0)
#' \dontrun{plot(fit1)}
#' 
#' # Mantel cross-correlogram 
#' fit2 <- mantel.correlog(dmat = dmat, zmat = zmat, wmat = wmat, increment = 2, 
#' quiet = TRUE, resamp = 0)
#' \dontrun{plot(fit2)}
#' @keywords spatial
#' @export
##############################################################################################
mantel.correlog<-function(dmat, zmat, wmat=NULL, increment, resamp = 1000, quiet=FALSE){
  ##############################################################################################
  
  if(is.null(wmat)){
    moran <- zmat[lower.tri(zmat)]
  }
  
  else {
    moran <- (zmat-mean(zmat))*(wmat-mean(wmat))/(sd(zmat)*sd(wmat))
    zero <- mean(diag(moran), na.rm= TRUE)
    moran <- moran[row(moran)>col(moran)]
  }
  
  if(resamp != 0){
    dmat2 <- dmat
    moran2 <- moran
  }
  
  if(is.null(wmat)){
    dmat <- dmat[lower.tri(dmat)]
  }
  else {
    dmat <- dmat[row(dmat)>col(dmat)]
  }
  
  
  dkl <- ceiling(dmat/increment)
  nlok <- sapply(split(moran, dkl), length)
  dmean <- sapply(split(dmat, dkl), mean, na.rm = TRUE)
  moran <- sapply(split(moran, dkl), mean, na.rm = TRUE)
  
  ly <- 1:length(dmean)
  x <- c(dmean[ly[moran < 0][1]], dmean[ly[moran < 0][1] - 1])
  y <- c(moran[ly[moran < 0][1] - 1], moran[ly[moran < 0][1]])
  if(moran[1] < 0) {
    tmp <- 0
  }
  else {
    if(sum(moran<0)>0){
      tmp <- lm(x ~ y)[[1]][1]
    }
    else{
      tmp <- NA
    }
  }
  
  p<-NULL
  
  if(resamp != 0){
    perm <- matrix(NA, ncol = length(moran), nrow = resamp)
    
    n<-dim(zmat)[1]
    for(i in 1:resamp){
      trekk <-sample(1:n)
      dma <- dmat2[trekk,trekk]
      mor <- moran2
      
      if(is.null(wmat)){
        dma <- dma[lower.tri(dma)]
      }
      else {
        dma <- dma[row(dma)>col(dma)]
      }
      
      dkl <- ceiling(dma/increment)	#generates the distance matrices
      perm[i,] <- sapply(split(mor, dkl), mean, na.rm = TRUE)
      
      whn=pretty(c(1,resamp), n=10)
      if(! quiet & any(i==whn)){
        cat(i, " of ", resamp, "\r")
        flush.console()}
    }
    
    p=(apply(moran<=t(perm),1,sum))/(resamp+1)
    p=apply(cbind(p, 1-p), 1, min) + 1/(resamp+1)
  }
  
  res <- list(n = nlok, mean.of.class = dmean, correlation = moran,
              x.intercept = tmp, p=p, call=deparse(match.call()))
  
  if(!is.null(wmat)){
    res$corr0 <- zero
  }
  class(res) <- "correlog"
  res
}


