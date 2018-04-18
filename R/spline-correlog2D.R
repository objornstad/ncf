#' @title Anisotropic nonparametric (cross-)correlation function for univariate spatial data
#' @description \code{spline.correlog2D} is the function to estimate the anisotropic nonparametric correlation function in 8 (or arbitrary) directions (North - Southeast) for univariate data. Correlation functions are calculated for each different bearing. The function assumes univariate observations at each location. (use \code{\link{Sncf2D}} otherwise).
#' @param x vector of length n representing the x coordinates.
#' @param y vector of length n representing the y coordinates.
#' @param z vector of length n representing the observation at each location.
#' @param w an optional second vector of length n for variable 2 (to estimate spatial or lagged cross-correlation functions).
#' @param df degrees of freedom for the spline. Default is sqrt(n).
#' @param type takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param npoints the number of points at which to save the value for the spline function (and confidence envelope / null distribution).
#' @param save If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).
#' @param max.it the maximum iteration for the Newton method used to estimate the intercepts.
#' @param xmax If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values for each pair of time series -- it will dump if any one pair has less than two (temporally) overlapping observations.
#' @param jitter If TRUE, jitters the distance matrix to avoid problems associated with fitting the function to data on regular grids.
#' @param quiet If TRUE, the counter is supressed during execution.
#' @param angle specifies number of cardinal directions and angles for which to calculate correlation functions. Default are 8 directions between 0 and 180.
#' @return An object of class "Sncf2D" is returned. See \code{\link{Sncf2D}} for details.
#' @details see \code{\link{Sncf2D}}
#' @note The function to estimate the UNIvariate anisotropic nonparametric (cross-)correlation function in arbitrary directions. In particular it was developed to calculate the univariate lagged cross-correlation function used in (Humston et al. 2005). Note that this 2D spline correlogram does the anisotopic analysis NOT by doing the angle-with-tolerance-wedge-style of Oden and Sokal (1986) but by projecting the the spatial coordinates of all locations on a sequence of cardinal angles (a la Sncf2D). Hence, all data points are used every time, it is only their relative distances that are changed. For example \{0, 0\} and \{0, 10\} are distance zero in the zero-degree direction but at distance 10 in the 90-degree direction.
#' @references Oden, N.L. and Sokal, R.R. 1986. Directional autocorrelation: an extension of spatial correlograms to two dimensions. Systematic Zoology 35: 608-617. \url{https://doi.org/10.2307/2413120}
#' 
#'   Humston, R., Mortensen, D. and Bjornstad, O.N. 2005. Anthropogenic forcing on the spatial dynamics of an agricultural weed: the case of the common sunflower. Journal of Applied Ecology 42: 863-872. \url{https://doi.org/10.1111/j.1365-2664.2005.01066.x}
#' @seealso \code{\link{Sncf2D}}
#' @keywords smooth regression
#' @export
###############################################################################################
spline.correlog2D <- function(x, y, z, w=NULL, df = NULL, type = "boot", resamp = 1000, npoints = 300,
                               save = FALSE, max.it=25, xmax=FALSE, na.rm = FALSE, jitter=FALSE, quiet=FALSE,
                               angle=c(0,22.5,45,67.5,90,112.5,135,157.5)){
  ##############################################################################################
  #spline.correlog2D is the function to estimate the anisotropic nonparametric covariance function
  #(using a smoothing spline as an equivalent kernel) in 8 (or arbitrary) directions (North - Southeast) 
  #through calculateing projected distances onto the different bearings (i.e. all data are used for each 
  #direction = 0,22.5,45,67.5,90,112.5,135,157.5)
  ############################################################################################
  
  #the following sets up the output:
  
  real<-lapply(unlist(lapply(angle, as.name)),as.null)
  names(real)<-unlist(lapply(angle, as.name))
  
  for(i in 1:length(angle)){
    real[[i]]<-list(x.intercept = NA, e.intercept = NA, y.intercept = NA, 
                    cbar.intercept = NA, predicted = list(x = matrix(NA, nrow = 1, 
                                                                     ncol = npoints), y = matrix(NA, nrow = 1, ncol = npoints)))
  }
  
  real$cbar<-NA
  
  if(resamp==0){
    boot<-lapply(unlist(lapply(angle, as.name)),as.null)
    names(boot)<-unlist(lapply(angle, as.name))
  }
  
  else{
    boot<-lapply(unlist(lapply(angle, as.name)),as.null)
    names(boot)<-unlist(lapply(angle, as.name))
    
    for(i in 1:length(angle)){
      boot[[i]]<-list(boot.summary=list(
        cbar = matrix(NA, ncol=1, nrow=resamp),
        x.intercept = matrix(NA, ncol=1, nrow=resamp),
        e.intercept = matrix(NA, ncol=1, nrow=resamp),
        y.intercept = matrix(NA, ncol=1, nrow=resamp),
        cbar.intercept = matrix(NA, ncol=1, nrow=resamp)
      ),
      predicted= list(x = matrix(NA, nrow = 1, ncol = npoints),
                      y = matrix(NA, nrow = resamp, ncol = npoints)))
    }
  }
  
  type <- charmatch(type, c("boot", "perm"), 
                    nomatch = NA)
  if(is.na(type))
    stop("method should be \"boot\", or \"perm\"")
  
  NAO <- FALSE
  
  #check for missing values
  if(any(!is.finite(z))) {
    if(na.rm){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pariwise deletion")
    }
  }
  
  #This generates the moran distances for cross-correlation
  #the odd adding of zero is just to ensure that all vectors 
  #are treated as numeric
  n <- length(z)
  #p <- dim(z)[2]
  #z <- as.matrix(z)+0
  
  #moran <- cor2(t(z), circ=FALSE)
  
  zscal <- (scale(z, center = TRUE, scale = TRUE)[, 1])/(sqrt((n-1)/n))
  
  if(is.null(w)){
    moran <- t(outer(zscal, zscal))
  }
  else {
    wscal <- (scale(w, center = TRUE, scale = TRUE)[, 1])/(sqrt((n-1)/n))
    zw <- c(zscal,wscal)
    moran <- t(outer(zw, zw))[1:n,(n+1):(2*n)]
  }
  
  
  if(is.null(df)){
    df <- sqrt(n)
  }
  
  maxdist <- ifelse(!xmax,max(sqrt(outer(x,x, "-")^2+outer(y,y,"-")^2)),xmax)
  xpoints <- seq(-maxdist, maxdist, length = npoints)
  
  #loop over directions
  ang <- (2*pi)*angle/360
  for(d in 1:length(ang)){
    #The next fits the spline function
    #then generating geographic distances
    y2 <- x * sin (ang[d]) + y * cos (ang[d])
    
    xdist <- outer(y2,y2,"-")
    
    if(jitter == TRUE) {
      #xdist <- jitter(xdist)
      xdist<-apply(xdist,2,jitter)
    }
    
    mdist <- max(xdist)
    
    if(is.null(w)){
      triang <- col(xdist)!=row(xdist)
    }
    
    else {
      triang <- xdist
      triang[,] <- TRUE
      triang <- triang==1
    }
    
    u <- xdist[triang]
    v <- moran[triang]
    sel <- is.finite(v) & is.finite(u)
    u <- u[sel]
    v <- v[sel]
    v <- v[abs(u) <= maxdist]
    u <- u[abs(u) <= maxdist]
    
    out=gather(u=u,v=v,w=w, moran=moran, df=df, xpoints=xpoints, filter=FALSE, fw=0)
    real$cbar <- mean(v)
    if(is.null(w)){
      real[[d]]$y.intercept <- out$yint
    }
    else{
      real[[d]]$y.intercept <-  mean(diag(moran))
    }
    
    real[[d]]$predicted <- list(x = out$x, y = out$y)
    real[[d]]$predicted$y[abs(real[[d]]$predicted$x)>mdist] <- NA
    real[[d]]$x.intercept <-out$xint
    real[[d]]$e.intercept <- out$eint
    real[[d]]$cbar.intercept <- out$cint
    
    ##End of spline fit
    
    if(resamp != 0) {
      #here is the bootstrapping/randomization
      boot[[d]]$predicted$x[1,] <- xpoints
      
      for(i in 1:resamp) {
        whn=pretty(c(1,resamp), n=10)
        if(! quiet & any(i==whn)){
          cat(i, " of ", resamp, "(direction", d, "of ", length(ang),")\r")
          flush.console()}
        if(type == 1) {
          trekkx <- sample(1:n, replace = TRUE)
          trekky <- trekkx
        }
        
        if(type == 2) {
          trekky <- sample(1:n, replace = FALSE)
          trekkx <- 1:n
        }
        
        xdistb <- xdist[trekkx, trekkx]
        
        if(is.null(w)){
          triang <- col(xdistb)!=row(xdistb)
        }
        
        else {
          triang <- xdistb
          triang[,] <- TRUE
          triang <- triang==1
        }
        
        xdistb <- xdistb[triang]
        moranb <- moran[trekky, trekky][triang]
        
        if(type == 1&is.null(w)) {
          moranb <- moranb[!(xdistb == 0)]
          xdistb <- xdistb[!(xdistb == 0)]
        }
        
        u <- xdistb
        v <- moranb
        sel <- is.finite(v) & is.finite(u)
        u <- u[sel]
        v <- v[sel]
        v <- v[u<=maxdist]
        u <- u[u<=maxdist]
        
        out=gather(u=u,v=v,w=w, moran=moranb, df=df, xpoints=xpoints, filter=FALSE, fw=0)
        boot[[d]]$boot.summary$cbar[i,1] <- mean(v)		
        
        if(is.null(w)){
          boot[[d]]$boot.summary$y.intercept[i,1] <- out$yint
        }
        else {
          boot[[d]]$boot.summary$y.intercept[i,1] <- mean(diag(moran[trekky, trekky]))
        }
        
        boot[[d]]$predicted$y[i,] <- out$y
        boot[[d]]$predicted$y[i,][abs(boot[[d]]$predicted$x[1,])>mdist] <- NA
        boot[[d]]$boot.summary$x.intercept[i,1] <- out$xint
        boot[[d]]$boot.summary$e.intercept[i,1] <- out$eint
        boot[[d]]$boot.summary$cbar.intercept[i,1]  <- out$cint
      }
      
      ##end of bootstrap loop!
      
      if(save == TRUE) {
        boot[[d]]$boot <- list(predicted = boot[[d]]$predicted)
      }
      else {
        boot[[d]]$boot <- NULL
      }
      
      ty <- apply(boot[[d]]$predicted$y, 2, quantile, probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 
                                                                0.75, 0.9, 0.95, 0.975, 1), na.rm = TRUE)
      dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), NULL)
      tx <- boot[[d]]$predicted$x
      boot[[d]]$boot.summary$predicted <- list(x = tx,y = ty)
    }
    
  }	
  res <- list(real = real, boot = boot, max.distance = maxdist, angle=angle, call=deparse(match.call()))
  class(res) <- "Sncf2D"
  res
}