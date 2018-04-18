#' @title The mean (cross-)correlation (with bootstrapp CI) for a panel of spatiotemporal data
#' @description \code{mSynch} is the function to estimate the mean (cross-)correlation in a spatiotemporal dataset as discussed in Bjornstad et al. (1999). The function requires multiple observations at each location.
#' @param x matrix of dimension n x p representing p observation at each location (i.e. each row is a time series).
#' @param y optional matrix of dimension m x p representing p observation at each location (i.e. each row is a time series). If provided, the mean cross-correlation between the two panels is computed.
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values for each pair of time series -- it will dump if any one pair has less than two (temporally) overlapping observations.
#' @param circ If TRUE, the observations are assumed to be angular (in radians), and circular correlation is used.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @return An object of class "mSynch" is returned, consisting of a list with two components: 
#' \item{real}{the regional average correlation.}
#' \item{boot}{a vector of bootstrap resamples.} 
#' @details Missing values are allowed -- values are assumed missing at random. 
#' 
#'   The circ argument computes a circular version of the Pearson's product moment correlation (see \code{\link{cor2}}).
#' @references Bjornstad, O.N., Ims, R.A. & Lambin, X. (1999) Spatial population dynamics: Analysing patterns and processes of population synchrony. Trends in Ecology and Evolution, 11, 427-431. \url{https://doi.org/10.1016/S0169-5347(99)01677-8}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{print.mSynch}}
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
#' # mean correlation analysis
#' fit1 <- mSynch(x = z, resamp = 500)
#' print(fit1)
#' @keywords spatial
#' @export
##############################################################################################
mSynch<-function(x, y=NULL, resamp = 1000, na.rm = FALSE, circ=FALSE, quiet=FALSE){
  ##############################################################################################
  #mSynch is a function to estimate the mean (cross-)correlation with bootstrapp CI for one
  #or two panels of spatiotemporal data
  ############################################################################################
  
  NAO <- FALSE
  
  #check for missing values
  if(any(!is.finite(x))) {
    if(na.rm){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  
  #the following sets up the output:
  Sbar <- list(real = NA)
  
  #first generates the correlations
  #the odd adding of zero is just to ensure that all vectors 
  #are treated as numeric
  n <- dim(x)[1]
  p <- dim(x)[2]
  x <- as.matrix(x)+0
  
  if(!is.null(y)){
    m <- dim(y)[1]
    y <- as.matrix(y)+0
  }
  
  if(!is.null(y)){
    synch <- cor2(t(x), t(y), circ=circ)
  }
  else{
    synch <- cor2(t(x), circ=circ)
  }
  
  if(is.null(y)){
    triang <- lower.tri(synch)
    v <- synch[triang]
  }
  else{
    v <- synch
  }
  
  Sbar$real <- mean(v, na.rm= TRUE)
  
  #End of real fit
  
  Sbar$boot <- NULL
  
  if(resamp != 0) {
    Sbar$boot<-matrix(NA, nrow = resamp, ncol = 1)
    for(i in 1:resamp) {
      whn=pretty(c(1,resamp), n=10)
      if(! quiet & any(i==whn)){
        cat(i, " of ", resamp, "\r")
        flush.console()}
      
      #here is the bootstrapping/randomization
      trekkx <- sample(1:n, replace = TRUE)
      
      if(!is.null(y)){
        trekky <- sample(1:m, replace = TRUE)
        synchb <- synch[trekkx,trekky]
      }
      
      else{
        trekky <- trekkx
        synchb <- synch[trekkx,trekkx][diag(n)[trekkx,trekkx]!=1]
      }
      
      Sbar$boot[i,1] <- mean(synchb, na.rm= TRUE)
    }
    #end of bootstrap loop!
  }
  res <- Sbar
  class(res) <- "mSynch"
  res
}

#' @title Print function for mSynch objects
#' @description `print' method for class "mSynch".
#' @param x an object of class "mSynch", usually, as a result of a call to \code{mSynch}.
#' @param verbose If TRUE, a raw listing of the object is produced. If FALSE, a summary list is produced
#' @param \dots other arguments
#' @return If verbose is FALSE, a list summarizing the regional correlation is produced: 
#' \item{mean}{the regional mean correlation.} 
#' \item{Squantile}{the quantile distribution from the resampling for the regional correlation.}
#' @seealso \code{\link{mSynch}}
#' @keywords smooth regression
##############################################################################################
"print.mSynch" <- function(x, verbose = FALSE, ...){
  ##############################################################################################
  #this is the generic print function for mSynch objects
  #
  #ARGUMENTS
  #verbose   If TRUE, the raw list is echoed
  ##############################################################################################
  if(!verbose) {
    Sbard <- quantile(x$boot, probs = c(0, 0.025, 0.25, 0.5,
                                        0.75, 0.975, 1), na.rm= TRUE)
    out <- list(mean=x$real, Squantile = Sbard)
    print(out)
    cat("\n\nFor a raw listing use print(x, verbose= TRUE)\n")
  }
  if(verbose) {
    print.default(x)
  }
}
