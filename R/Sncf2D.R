#' @title Anisotropic nonparametric (cross-)correlation function for spatio-temporal data
#' @description \code{Sncf2D} is the function to estimate the anisotropic nonparametric correlation function in 8 (or arbitrary) directions (North - Southeast). Correlation functions are calculated for each different bearing. The function requires multiple observations at each location. (use \code{\link{spline.correlog.2D}} otherwise).
#' @param x vector of length n representing the x coordinates.
#' @param y vector of length n representing the y coordinates.
#' @param z matrix of dimension n x p representing p observation at each location.
#' @param w an optional second matrix of dimension n x p for variable 2 (to estimate spatial or lagged cross-correlation functions).
#' @param df degrees of freedom for the spline. Default is sqrt(n).
#' @param type takes the value "boot" (default) to generate a bootstrap distribution or "perm" to generate a null distribution for the estimator
#' @param resamp the number of resamples for the bootstrap or the null distribution.
#' @param npoints the number of points at which to save the value for the spline function (and confidence envelope / null distribution).
#' @param save If TRUE, the whole matrix of output from the resampling is saved (an resamp x npoints dimensional matrix).
#' @param max.it the maximum iteration for the Newton method used to estimate the intercepts.
#' @param xmax If FALSE, the max observed in the data is used. Otherwise all distances greater than xmax is omitted.
#' @param na.rm If TRUE, NA's will be dealt with through pairwise deletion of missing values for each pair of time series -- it will dump if any one pair has less than two (temporally) overlapping observations.
#' @param jitter If TRUE, jitters the distance matrix, to avoid problems associated with fitting the function to data on regular grids
#' @param quiet If TRUE, the counter is supressed during execution.
#' @param angle specifies number of cardinal directions and angles for which to calculate correlation functions. Default are 8 directions between 0 and 180.
#' @return An object of class "Sncf2D" is returned, consisting of a list of estimates for each cardinal direction :
#' \item{real}{the list of estimates from the data.}
#' \item{$cbar}{the regional average correlation.}
#' \item{$x.intercept}{the lowest value at which the function is = 0. If correlation is initially negative, the distance is given as negative.}
#' \item{$e.intercept}{the lowest value at which the function 1/e.}
#' \item{$y.intercept}{the extrapolated value at x=0 (nugget).}
#' \item{$cbar.intercept}{distance at which regional average correlation is reach.}
#' \item{$predicted$x}{the x-axes for the fitted covariance function.}
#' \item{$predcited$y}{the values for the covariance function.}
#' \item{boot}{a list with the analogous output from the bootstrap or null distribution.} 
#' \item{$summary}{gives the full vector of output for the x.intercept, y.intercept, e.intercept, cbar.intercept, and the cbar and a quantile summary for the resampling distribution.}
#' \item{$boot}{If save=TRUE, the full raw matrices from the resampling is saved.}
#' \item{angle}{a vector with the cardinal directions.}
#' \item{max.distance}{the maximum spatial distance.}
#' @details Correlation functions are calculated on projected distances onto the different bearings so ALL data are used for each direction. The (obsolete?) \code{oldncf2D} used the alternative of slizing up the data like pieces of a pie.
#' 
#'   Latitude-longitude coordinates can NOT be used.
#'   
#'   Missing values are allowed - values are assumed missing at random.
#'   
#'   I have implemented an optional argument: \code{jitter} if TRUE this jitters the distance matrix, to avoid some problems I've had with spline-smoothing data from regular grid-data.
#' @note The function to estimate the anisotropic nonparametric (cross-)correlation function in arbitrary directions. In particular it was developed to calculate the lagged cross-correlation function (Bjornstad et al. 2002).
#' @references Bjornstad, O. N., M. Peltonen, A. M. Liebhold, and W. Baltensweiler. 2002. Waves of larch budmoth outbreaks in the European Alps. Science 298:1020-1023. \url{https://doi.org/10.1126/science.1075182}
#' @author Ottar N. Bjornstad \email{onb1@psu.edu}
#' @seealso \code{\link{summary.Sncf2D}}, \code{\link{plot.Sncf2D}}, \code{\link{cc.offset}} , \code{\link{Sncf}}, \code{\link{spline.correlog.2D}}
#' @examples 
#' # first generate some sample data
#' x <- expand.grid(1:20, 1:5)[, 1]
#' y <- expand.grid(1:20, 1:5)[, 2]
#' # z data from an exponential random field
#' z <- cbind(
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp"),
#'   rmvn.spa(x = x, y = y, p = 2, method = "exp")
#'   )
#' # anisotorpic nonparametric covariance function at 30 and 60 degrees
#' fit1 <- Sncf2D(x = x, y = y, z = z, resamp = 0, angle = c(30, 60))
#' \dontrun{plot.Sncf2D(fit1)}
#' summary.Sncf2D(fit1)
#' 
#' # What distance is the peak in correlation
#' cc.offset(fit1)
#' @keywords smooth regression
##############################################################################################
Sncf2D <- function(x, y, z, w=NULL, df = NULL, type = "boot", resamp = 1000, npoints = 300,
                   save = FALSE, max.it=25, xmax= FALSE, na.rm = FALSE, jitter= FALSE, quiet=FALSE,
                   angle=c(0,22.5,45,67.5,90,112.5,135,157.5)){
  ##############################################################################################
  #Sncf2D is the function to estimate the anisotropic nonparametric covariance function 
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
        cbar.intercept = matrix(NA, ncol=1, nrow=resamp)),
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
  if(any(!is.finite(unlist(z)))) {
    if(na.rm){
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pariwise deletion")
    }
  }
  
  if(is.null(w)){
    #This generates the moran distances for cross-correlation
    #the odd adding of zero is just to ensure that all vectors 
    #are treated as numeric
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
      #this is corrected to be
      #xdist <- jitter(xdist)
      xdist <- apply(xdist,2,jitter)
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
        if(! quiet & any(i==whn)){cat(i, " of ", resamp, "(direction", d, "of ", length(ang),")\r")
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

#' @title Print function for Sncf2D objects
#' @description `print' method for class "Sncf2D".
#' @param x an object of class "Sncf2D", usually, as a result of a call to \code{Sncf2D} or \code{spline.correlog.2D}).
#' @param \dots other arguments
#' @return The function-call is printed to screen.
#' @seealso \code{\link{Sncf2D}}
##############################################################################################
print.Sncf2D=function(x, ...){
  ##############################################################################################
  cat("This is an object of class Sncf2D produced by the call:\n\n", x$call, "\n\n Use summary() or plot() for inspection,  (or print.default() to see all the gory details).")}

#' @title Plots anisotropic spatial correlation-functions
#' @description `plot' method for class "Sncf2D".
#' @param x an object of class "Sncf2D", ususally, as a result of a call to \code{Sncf2D}.
#' @param xmax the maximal distance to be plotted on the x-axis. If set to zero the maximum distance in the data will be used.
#' @param ylim limits for the y-axis (default: -1, 1).
#' @param detail If TRUE, a separate plot is made for each direction (including confidence envelopes; see \code{\link{plot.Sncf}} for details. If FALSE, all correlation functions are superimposed on the same plot.
#' @param \dots other arguments
#' @return A plot or panel-plot results. These represents the xy-plot of distance against spatial (cross-)correlation for each cardinal direction.
#' @seealso \code{\link{Sncf2D}}, \code{\link{plot.Sncf}}
#' @keywords smooth regression
##############################################################################################
plot.Sncf2D <- function(x, xmax = 0, ylim=c(-1,1), detail = FALSE, ...){
  ##############################################################################################
  #this is the generic plot function for Sncf2D objects
  ##############################################################################################
  L<-length(x$angle)
  
  xmax <- ifelse(xmax == 0, x$max.distance, xmax)
  plot(x$real[[1]]$predict$x, x$real[[1]]$predict$y, xlim = c(-xmax, xmax), ylim
       = ylim, type = "n", xlab = "", ylab = "")
  lines(c(-max(x$real[[1]]$predict$x), max(x$real[[1]]$predict$x)), c(0, 0))
  lines(c(-max(x$real[[1]]$predict$x), max(x$real[[1]]$predict$x)), c(x$real$cbar, x$real$cbar))
  
  if(detail){
    
    par(mfrow=c(ceiling(sqrt(L)),ceiling(sqrt(L))))
    
    for(i in 1:L){
      plot(x$real[[i]]$predict$x, x$real[[i]]$predict$y, xlim = c(-xmax, xmax), ylim
           = ylim, type = "l", xlab = "Distance", ylab = "Correlation")
      
      if(!is.null(x$boot[[i]]$boot.summary)){
        polygon(c(x$boot[[i]]$boot.summary$predicted$x,rev(x$boot[[i]]$boot.summary$predicted$x)), c(x$boot[[i]]$boot.summary$predicted$y["0.025",], rev(x$boot[[i]]$boot.summary$predicted$y["0.975",])), col=gray(0.8), lty=0)
      }	
      lines(x$real[[i]]$predict$x, x$real[[i]]$predict$y)
      lines(c(-max(x$real[[i]]$predict$x), max(x$real[[i]]$predict$x)), c(0, 0))
      lines(c(-max(x$real[[i]]$predict$x), max(x$real[[i]]$predict$x)), c(x$real$cbar, x$real$cbar))
    }
  }
}

#' @title Summarizing anisotropic spatial correlation-functions
#' @description `summary' method for class "Sncf2D".
#' @param object an object of class "Sncf2D", usually, as a result of a call to \code{\link{Sncf2D}}.
#' @param \dots other arguments
#' @return A list summarizing the nonparametric covariance function in each cardinal direction results, each with the entires as in \code{\link{summary.Sncf}}.
#' @seealso \code{\link{Sncf2D}}, \code{\link{cc.offset}}, \code{\link{summary.Sncf}}
#' @keywords smooth regression
##############################################################################################
summary.Sncf2D<-function(object, ...){
  ##############################################################################################
  #this is the generic summary function for Sncf objects
  ##############################################################################################
  L<-length(object$angle)
  
  xy <- matrix(NA, ncol=L, nrow=4)
  dimnames(xy) <- list(c("x", "e", "cbar", "y"),
                       unlist(lapply(object$angle, as.name)))
  xyd <- lapply(unlist(lapply(object$angle, as.name)),as.null)
  names(xyd)<-unlist(lapply(object$angle, as.name))
  
  for(i in 1:L){
    xyd[[i]]<-matrix(NA, ncol=7, nrow=4)
  }
  
  cbar <- round(object$real$cbar, 2)
  
  for(i in 1:L){
    xy[1,i] <- round(object$real[[i]]$x.intercept, 2)
    xy[2,i] <- round(object$real[[i]]$e.intercept, 2)
    xy[3,i] <- round(object$real[[i]]$cbar.intercept, 2)
    xy[4,i] <- round(object$real[[i]]$y.intercept, 2)
    if(!is.null(object$boot[[i]]$boot.summary)){
      yd <- apply(object$boot[[i]]$boot.summary$y.intercept, 2, quantile, probs = c(0, 0.025, 0.25, 0.5,
                                                                                    0.75, 0.975, 1), na.rm= TRUE)
      xd <- apply(object$boot[[i]]$boot.summary$x.intercept, 2, quantile, probs = c(0, 0.025, 0.25, 0.5,
                                                                                    0.75, 0.975, 1), na.rm= TRUE)
      ed <- apply(object$boot[[i]]$boot.summary$e.intercept, 2, quantile, probs = c(0, 0.025, 0.25, 0.5,
                                                                                    0.75, 0.975, 1), na.rm= TRUE)
      synchd <- quantile(object$boot[[i]]$boot.summary$cbar[,1], probs = c(0, 0.025, 0.25, 0.5,
                                                                           0.75, 0.975, 1), na.rm= TRUE)
      cbard <- quantile(object$boot[[i]]$boot.summary$cbar.intercept[,1], probs = c(0, 0.025, 0.25, 0.5,
                                                                                    0.75, 0.975, 1), na.rm= TRUE)
      xyd[[i]] <- cbind(xd, ed, yd, cbard)
      dimnames(xyd[[i]]) <- list(c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), c("x", "e", "y", "cbar"))}
    if(is.null(object$boot[[i]]$boot.summary)){
      synchd <- NULL
      xyd <- NULL	
    }
  }
  res <- list(call = object$call, Regional.synch = object$real$cbar, Squantile = synchd, estimates = xy, quantiles = xyd)
  res
}

#' @title Function to calculate the distance at which the cross-correlation peaks for Sncf objects
#' @description Alternative `summary' method for class "Sncf2D".
#' @param object an object of class "Sncf2D", usually, as a result of a call to \code{Sncf2D} or \code{spline.correlog.2D}.
#' @param xmax the maximum distance to consider (default is no upper limit).
#' @return An matrix of class "cc.offset" is returned with columns:
#' \item{angle}{the cardinal angle (in degrees).}
#' \item{distance}{the distances (in the positive direction) to the mode of the (cross-) correlation function (with 95\% confidence bounds).}
#' \item{correlation}{the correlation at the mode (with CI) for each of the cardinal angles.}
#' @seealso \code{\link{Sncf2D}}, \code{\link{summary.Sncf2D}}, \code{\link{plot.cc.offset}}
#' @keywords smooth regression
##############################################################################################
cc.offset<-function(object, xmax=NULL){
  ##############################################################################################
  #calculates offsets in Sncf2D cross-correltions
  ##############################################################################################
  
  L<-length(object$angle)
  
  if(is.null(xmax)){
    xmax<-object$max.distance
  }
  
  if(is.null(object$boot[[1]]$boot.summary)){
    xy <- matrix(NA, ncol=3, nrow=L)
    dimnames(xy) <- list(unlist(lapply(object$angle, as.name)),
                         c("angle", "distance", "correlation"))
  }
  if(!is.null(object$boot[[1]]$boot.summary)){
    xy <- matrix(NA, ncol=7, nrow=L)
    dimnames(xy) <- list(unlist(lapply(object$angle, as.name)),
                         c("angle", "distance", "correlation", "dL", "dU", "cL", "cU"))
  }
  
  xy[,1]<-object$angle
  
  for(i in 1:L){
    sel<-object$real[[i]]$predict$x>=0&object$real[[i]]$predict$x<xmax
    D<-object$real[[i]]$predict$x[sel]
    D2<-object$real[[i]]$predict$y[sel]
    xy[i,2]<-na.omit(D[D2==max(D2, na.rm= TRUE)])[1]
    xy[i,3]<-na.omit(D2[D2==max(D2, na.rm= TRUE)])[1]
    if(!is.null(object$boot[[i]]$boot.summary)){
      xy[i,4:5]<-range(object$boot[[i]]$boot.summary$predict$x[sel][(object$boot[[i]]$boot.summary$predict$y["0.975", sel]-max(D2, na.rm= TRUE))>0], na.rm= TRUE)
      xy[i, 7]<-object$boot[[i]]$boot.summary$predict$y["0.975", sel][rev(order(na.omit(D2)))[1]]
      xy[i, 6]<-object$boot[[i]]$boot.summary$predict$y["0.025", sel][rev(order(na.omit(D2)))[1]]
    }
  }
  class(xy)<-"cc.offset"
  return(xy)
}

#' @title Plots the cc.offset summary of the anisotropic spatial correlation-functions
#' @description `plot' method for class "cc.offset".
#' @param x an object of class "cc.offset", ususally, as a result of applying \code{cc.offset} to an object of class \code{Sncf2D}.
#' @param dmax the maximal distance for radial plot. If NULL, the maximum distance in the data will be used.
#' @param inches the size of the symbols.If NULL, default is 0.1.
#' @param \dots other arguments
#' @return A radial `symbol' plot results. The radius represents the distance to peak correlation (the mode) of the correlation function (in the positive direction). The size of the symbol represents the magnitude of the correlation at the mode for the given cardinal direction.
#' @seealso \code{\link{cc.offset}}, \code{\link{Sncf2D}}, \code{\link{plot.Sncf2D}}
##############################################################################################
plot.cc.offset <- function(x, dmax=NULL, inches=NULL, ...){
  ##############################################################################################
  #this is the generic plot function for cc.offset objects
  ##############################################################################################
  theta <- 2*pi*x[,"angle"]/360
  x2 <- x[,"distance"]*sin(theta)
  y <- x[,"distance"]*cos(theta)
  inc=inches
  if(is.null(dmax)) dmax=max(abs(c(x2,y)))
  if(is.null(inches)) inc=.1
  #tmp<-rep(0, length(theta))
  tmp<-x[,"correlation"]
  axs=pretty(c(0, dmax), n=4)
  yl=xl=c(-max(axs), max(axs))
  symbols(x2,y, circles= ifelse(tmp>0,tmp,0), inches=inc, xlim=xl, ylim=yl, xlab="", ylab="", fg=1, bg=2, asp=1, xaxt="n", yaxt="n", bty="n")
  symbols(rep(0,length(axs)), rep(0,length(axs)), circles=axs, inches=FALSE, xlab="", 
          ylab="", xaxt="n", yaxt="n", bty="n", add=TRUE)	
  lines(c(0,0), yl)
  lines(yl, c(0,0))
  text(axs[-1], rep(0,length(axs))[-1], axs[-1], cex=0.6, pos=1)
  symbols(x2,y, circles= ifelse(tmp>0,tmp,0), inches=inc, xlim=xl, ylim=yl, xlab="", ylab="", fg=1, bg=2, asp=1, xaxt="n", yaxt="n", add=TRUE, bty="n")
}