#' @title Fit a Sum of Smooth Exponentials
#'
#' @description This function estimates a three-component non-parametric model of the age-specific mortality rates.
#' Its design resembles the parametric Heligman-Pollard model but it uses \emph{P}-splines instead of parametric functions to describe the force of mortality.
#'
#'
#' @param data data frame produced with \code{HMD2MH} or similarly structured
#' @param lambda.hump smoothing parameter for the hump component
#' @param lambda.sen smoothing parameter for the senescence part
#' @param x.hump upper age to compute the starting values for the hump component
#' @param x.sen lower age to compute the starting values for the senescence component
#' @param kappa ridge penalty
#' @param maxit maximum number of iterations of the optimization algorithm
#' @param verbose should the function print information during the optimization
#'
#'
#' @details Three non-parametric components are estimated for each phase of the force of mortality.
#' A first component estimates the decreasing trend observed during the first years of life.
#' A second component estimates the concave pattern observed during adolescence and early adulthood.
#' A third component estimates the approximately exponential increase observed during late adulthood.
#'
#' The function first estimates separate \emph{P}-splines for each components. Then, it uses a re-weighted iterative least square algorithm to adjust each component so that they sum up to the observed death rates.
#' The algorithm is relatively sensitive to starting values (i.e. the independent fit of each component). In case of a failed convergence, it is thus advised to reconsider the following parameters.
#'
#' The \code{x.hump} and \code{x.sen} arguments are used to compute the starting values of each component. More specifically, they are used to determined the age range on which the young adult mortality hump and the senescence components are fitted.
#' The values of \code{x.hump} and \code{x.sen} define in a way the interval where the hump and the senescence components overlap.
#' In case of a very narrow hump, it is advised to reduce the value of \code{x.hump} and/or of \code{x.sen}, wherease in case of an especially wide hump it may be useful to consider larger values for \code{x.hump} and \code{x.sen}.
#' Inadequate values of \code{x.hump} may result in incoherent starting values for the cause-deleted SSE models and a lack of convergence.
#'
#' The \code{lambda.hump} and \code{lambda.sen} parameters control the amount of smoothing imposed on the hump and senescence component respectively. Typically, the former is smaller than the latter, since this is the more instable part of the force of mortality.
#' However, in some cases, especially when using abridged datasets, it may be useful to consider smaller values of \code{lambda.sen}.
#'
#' The \code{maxit} argument defines the maximum number of iterations used for the adjustment step. This step uses a Penalized Composite Link Model (PCLM) along with an iterative re-weighted least squares (IRWLS) optimization procedure.
#' The \code{maxit} argument will therefore determine the maximum number of iterations used in the IRWLS loop, but also the speed of convergence since the re-weighting defines updated solution as \deqn{new = old * (1 - it/maxit) + new * (it/maxit)}
#' The \code{maxit} argument is defined at 200 by default, but it can be increased in case of an absence of convergence, even if the algorithm stopped before reaching \code{maxit} number of iterations.
#'
#' The \code{kappa} parameter regulates the force of the ridge penalty, which enforces the shape constraints, i.e. that the hump component is concave and the senescence component is monotonically increasing.
#'
#' @return
#' \describe{
#'    \item{data}{The original data frame produced with \code{HMD2MH} or similarly structured. Used later for the plot and summary methods.}
#'    \item{type}{The type of model used: non-parametric. To use in the summary method.}
#'    \item{method}{The algorithm used to fit the model: IRWLS. To use in the summary method.}
#'    \item{it}{The number of iterations of the IRWLS algorithm used to fit the model.}
#'    \item{maxit}{The maximum number of iterations of the optimization algorithm set before the estimation.}
#'    \item{deg2}{The degrees of freedom used for the \emph{P}-spline of the senescence component.}
#'    \item{deg3}{The degrees of freedom used for the \emph{P}-spline of the hump component.}
#'    \item{coef}{The concanated fitted coefficients of the \emph{P}-splines.}
#'    \item{XX}{List containing the base spline functions for each of the three components.}
#'    \item{mhat}{List containing three vectors containing each the age-specific mortality contributions of the estimated components. The sum of \code{mhat1}, \code{mhat2} and \code{mhat3} gives the overall fitted force of mortality.}
#' }
#'
#'
#' @examples
#'
#' data("CHE2010m")
#'
#' fit <- sse.fit(data = CHE2010m)
#'
#' @references
#'
#' Camarda, C. G., Eilers, P. H. C., & Gampe, J. (2016). Sums of smooth exponentials to decompose complex series of counts. Statistical Modelling.
#'
#' @seealso
#' \link{morthump}
#'
#' @export
#'
#' @import MortalitySmooth

morthump <- function(data, model, method = "port", w = 1/data$m, start = NULL, maxit = 500, x1 = 30, x2 = 50, lambda.sen = 100, lambda.hump = 1){
  
  if(model == "hp"){
    
    fit <- hp.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hp"
    
  }
  
  if(model == "hpk"){
    
    fit <- hpk.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hpk"
    
  }
  
  if(model == "hps"){
    
    fit <- hps.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hps"
    
  }
  
  if(model == "sse"){
    
    fit <- sse.fit(data, maxit = maxit, x.hump = x1, x.sen = x2, lambda.sen = lambda.sen, lambda.hump = lambda.hump)
    attr(fit,"model") <- "sse"
    
  }
  
  class(fit) <- "morthump"
  
  return(fit)
}

sse.fit <- function(data, lambda.hump = 1, lambda.sen = 10, x.hump = 35, x.sen = 50, kappa = 10^5, maxit = 200, verbose = FALSE){
  
  x <- data$x
  n <- length(x)
  x <- 1:n
  y <- data$d
  e <- data$n
  lmx <- log(y/e)
  
  ## FITTING SECTION
  
  staINF <- 1
  endINF <- n
  staACC <- 1
  endACC <- n
  staAGI <- 1
  endAGI <- n
  
  ## let each component work in a particular age-range
  ## infant mortality
  p1 <- c(staINF, endINF)
  x1 <- x[x>=p1[1] & x<=p1[2]]
  w1 <- as.numeric(x%in%x1)
  y1 <- y[which(w1==1)]
  e1 <- e[which(w1==1)]
  ## aging mortality
  p2 <- c(staAGI, endAGI)
  x2 <- x[x>=p2[1] & x<=p2[2]]
  w2 <- as.numeric(x%in%x2)
  y2 <- y[which(w2==1)]
  e2 <- e[which(w2==1)]
  ## accident mortality
  p3 <- c(staACC, endACC)
  x3 <- x[x>=p3[1] & x<=p3[2]]
  w3 <- as.numeric(x%in%x3)
  y3 <- y[which(w3==1)]
  e3 <- e[which(w3==1)]
  
  ## Design matrix for the infant part
  B1 <- cbind(1,1/x1)
  nb1 <- ncol(B1)
  ## difference matrices
  D1 <- diff(diag(nb1), diff=2)
  tDD1 <- t(D1) %*% D1
  
  ## B-splines for the aging part
  deg2 <- floor(length(x)/5)
  xl2 <- min(x)
  xr2 <- max(x)
  xmax2 <- xr2 + 0.01 * (xr2 - xl2)
  xmin2 <- xl2 - 0.01 * (xr2 - xl2)
  
  print(x)
  print(paste0('xMin2 = ',xmin2))
  print(paste0('xMax2 = ',xmax2))
  print(paste0('deg2 = ',deg2))
  
  B2 <- MortSmooth_bbase(x, xmin2, xmax2, deg2, 3)
  nb2 <- ncol(B2)
  ## difference matrices
  D2 <- diff(diag(nb2), diff=2)
  tDD2 <- t(D2) %*% D2
  
  ## B-splines for the accident-hump
  deg3 <- floor(length(x3)/5)
  xl3 <- min(x3)
  xr3 <- max(x3)
  xmax3 <- xr3 + 0.01 * (xr3 - xl3)
  xmin3 <- xl3 - 0.01 * (xr3 - xl3)
  
  print(x3)
  print(paste0('xMin3 = ',xmin3))
  print(paste0('xMax3 = ',xmax3))
  print(paste0('deg3 = ',deg3))
  B3 <- MortSmooth_bbase(x3, xmin3, xmax3, deg3, 3)
  nb3 <- ncol(B3)
  ## difference matrices
  D3 <- diff(diag(nb3), diff=3)
  tDD3 <- t(D3) %*% D3
  
  ## complete model matrix as a list
  ## in their order and with dimensions
  XX <- list(X1=B1, X2=B2, X3=B3)
  nx <- length(XX)
  nc <- unlist(lapply(XX, ncol))
  ind2 <- cumsum(nc)
  ind1 <- c(1, ind2[1:(nx-1)]+1)
  ind <- NULL # loop to avoid !
  for(i in 1:nx){
    ind[[i]] <- ind1[i]:ind2[i]
  }
  ## indicators for the fitted values
  indF <- cbind(w1, w2, w3)
  
  ## fitting a Poisson-GLM for infant mortality
  ww1 <- rep(0,length(y1))
  ww1[x1<=10] <- 1
  y1a <- c(y)
  y1a[which(y1a<=0)] <- 0
  y1 <- y1a[which(w1==1)]
  fit1 <- glm(y1~I(1/x1), offset=log(e1),
              weights=ww1,
              family=poisson)
  y1.st0 <- exp(XX[[1]] %*% fit1$coef)*e1
  lhaz1.st0 <- log(y1.st0/e1)
  
  y1.st <- rep(0,n)
  y1.st[which(w1==1)] <- y1.st0
  lhaz1.st <- log(y1.st/e)
  
  
  ## fitting P-splines over all x for the aging component
  y2a <- c(y-y1.st)
  y2a[which(y2a<=0)] <- 0
  y2 <- y2a[which(w2==1)]
  ww2 <- rep(0,length(y2))
  ww2[x2>=x.sen] <- 1
  fit2 <- suppressWarnings(Mort1Dsmooth(x=x2, y=y2, offset=log(e2),
                                        w=ww2,
                                        method=3, lambda=1000,
                                        ndx = deg2, deg = 3, pord = 2))
  y2.st0 <- exp(XX[[2]] %*% fit2$coef)*e2
  lhaz2.st0 <- log(y2.st0/e2)
  y2.st <- rep(0,n)
  y2.st[which(w2==1)] <- y2.st0
  lhaz2.st <- log(y2.st/e)
  
  ## fitting P-splines over all x for the accident component
  # if y1+y2 > y => y3 < 0 => no fit / solution = only if possible
  if(all(y1.st + y2.st > y)){
    y3a <- c(y-y1.st-y2.st)
    y3a[which(y3a<=0)] <- 0
    y3 <- y3a[which(w3==1)]}
  ww3 <- rep(0, length(x3))
  ww3[x3>15 & x3<x.hump] <- 1
  fit3 <- suppressWarnings(Mort1Dsmooth(x=x3, y=y3, offset=log(e3),
                                        w=ww3,
                                        method=3, lambda=100000,
                                        ndx = deg3, deg = 3, pord = 3))
  y3.st0 <- exp(XX[[3]] %*% fit3$coef)*e3
  lhaz3.st0 <- log(y3.st0/e3)
  y3.st <- rep(0,n)
  y3.st[which(w3==1)] <- y3.st0
  lhaz3.st <- log(y3.st/e)
  
  ## starting values for the whole fuction
  y.st <- y1.st + y2.st + y3.st
  lhaz.st <- log(y.st/e)
  
  ## concatenating the starting values
  coef.st <- as.vector(c(fit1$coef,
                         fit2$coef,
                         fit3$coef))
  coef <- coef.st
  
  ## penalty term just for the second component
  P2 <- lambda.sen * tDD2
  P3 <- lambda.hump * tDD3
  P <- bdiag(list(diag(0,nc[1]), P2, P3))
  
  ## including monotonicy for infant part
  D1mon <- diff(diag(nb1), diff=1)
  w1mon <- rep(0, nrow(D1mon))
  W1mon <- diag(w1mon)
  ## including monotonicy for aging part
  D2mon <- diff(diag(nb2), diff=1)
  w2mon <- rep(0, nrow(D2mon))
  W2mon <- diag(w2mon)
  ## including concaveness for accident-hump
  D3con <- diff(diag(nb3), diff=2)
  w3con <- rep(0, nrow(D3con))
  W3con <- diag(w3con)
  
  #plot(coef)
  
  ## iteration
  for(it in 1:maxit){
    
    d <- 1 - it/maxit
    
    ## penalty for the shape constraints
    P2mon <- kappa * t(D2mon) %*% W2mon %*% D2mon
    P3con <- kappa * t(D3con) %*% W3con %*% D3con
    Psha <- bdiag(list(diag(rep(0,nc[1])),
                       P2mon,
                       P3con))
    
    eta <- numeric(nx*n)
    
    
    for(i in 1:nx){
      eta0 <- rep(0, n)
      eta0[which(indF[,i]==1)] <- XX[[i]] %*% coef[ind[[i]]]
      eta[1:n+(i-1)*n] <- eta0
    }
  

    gamma <- exp(eta)*c(indF)
    
    
    

    
    mu <- numeric(n)
    for(i in 1:nx){
      mu <- (e * gamma[1:n+(i-1)*n]) + mu
    }
    
    w <- mu

    
    U <- matrix(NA, n, sum(nc))
    for(i in 1:nx){
      u <- gamma[1:n+(i-1)*n]/mu * e
      XXi <- matrix(0, nrow=n, ncol=nc[i])
      XXi[which(indF[,i]==1), ] <- XX[[i]]
      U0 <- u * XXi
      U[,ind[[i]]] <- U0
    }
    
    
  
    tUWU <- t(U) %*% (w * U)
    tUWUpP <- tUWU + P + Psha
    r <- y - mu
    tUr <- t(U) %*% r


    
    coef.old <- coef
    coef <- solve(tUWUpP, tUr + tUWU %*% coef)
   
    coef <- d * coef.old + (1 - d) * coef
    

    ## update weights for shape constraints
    ## aging
    W2mon.old <- W2mon
    coef2 <- coef[ind[[2]]]
    W2mon <- diag(diff(coef2) <= 0)
    ## accident-hump
    W3con.old <- W3con
    coef3 <- coef[ind[[3]]]
    W3con <- diag(diff(coef3, diff=2) >= 0)
    ## convergence criterion for the concaveness
    conv.con <- all(W2mon.old==W2mon,
                    W3con.old==W3con)
   #DEBUGGGGGG
    #print(it)
    
    dif.coef <- max(abs(coef.old - coef))/max(abs(coef))
    
    if(dif.coef < 1e-04 & it > 4 & conv.con) break
    if(verbose == T){cat(it, dif.coef, conv.con, "\n")}
  }
  
  
  coef.hat <- coef
  ## fitted values
  etas <- NULL
  for(i in 1:nx){
    etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
  }
  eta1.hat <- etas[[1]]
  eta2.hat <- etas[[2]]
  eta3.hat <- etas[[3]]
  
  eta.hat <- numeric(nx*n)
  for(i in 1:nx){
    eta0 <- rep(0, n)
    eta0[which(indF[,i]==1)] <- XX[[i]]%*%coef.hat[ind[[i]]]
    eta.hat[1:n+(i-1)*n] <- eta0
  }
  
  gamma1.hat <- exp(eta1.hat)
  gamma2.hat <- exp(eta2.hat)
  gamma3.hat <- exp(eta3.hat)
  gamma.hat <- exp(eta.hat)*c(indF)
  
  mu1.hat <- exp(eta1.hat)*e1
  mu2.hat <- exp(eta2.hat)*e
  mu3.hat <- exp(eta3.hat)*e3
  mu.hat <- numeric(n)
  for(i in 1:nx){
    mu.hat <- (e * gamma.hat[1:n+(i-1)*n]) + mu.hat
  }
  
  fit <- list(data = data, type = "non-parametric", method = "IRWLS", it = it, maxit = maxit, deg2 = deg2, deg3 = deg3, coef = coef.hat, XX = XX, mhat = list(mhat1 = gamma1.hat, mhat2 = gamma2.hat, mhat3 = gamma3.hat))
  
  return(fit)
  
}


sse2.fit <- function(x, d, e, lambda1, lambda2, kappa = 10^6, deg = NULL,
                     plotIT = TRUE, plotFIT = TRUE, lab = " ", mon = TRUE,
                     max.it=200, ridge=10^-4,
                     x1 = 35, x2 = 50){
  
  n <- length(x)
  y <- d
  lmx <- log(y/e)
  
  wei <- w1 <- w2 <- rep(1,n)
  wei[e==0] <- 0
  lhaz.act <- log(y/e)
  
  ## B-splines stuff common to both component
  if(is.null(deg)){
    deg <- floor(length(x)/5)
  }else{
    deg=deg
  }
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  B <- MortSmooth_bbase(x, xmin, xmax, deg, 3)
  nb <- ncol(B)
  ## difference matrices for accident-hump
  D1 <- diff(diag(nb), diff=3)
  tDD1 <- t(D1) %*% D1
  ## difference matrices for accident-hump
  D2 <- diff(diag(nb), diff=2)
  tDD2 <- t(D2) %*% D2
  
  ## complete model matrix as a list
  ## in their order and with dimensions
  XX <- list(X1=B, X2=B)
  nx <- length(XX)
  nc <- unlist(lapply(XX, ncol))
  ## indicators for the coefficients of each basis
  ind2 <- cumsum(nc)
  ind1 <- c(1, ind2[1:(nx-1)]+1)
  ind <- NULL
  for(i in 1:nx){
    ind[[i]] <- ind1[i]:ind2[i]
  }
  ## indicators for the fitted values
  indF <- cbind(w1, w2)
  ## number of total coefficients
  np <- sum(nc)
  
  ## starting values
  
  ## fitting P-splines for the aging
  ## taking only ages "x.sen+" and extrapolating
  ## for the previous ages
  y2 <- y[which(w2==1)]
  ww2 <- rep(0,length(y2))
  ww2[x>=x2] <- 1
  fit2 <- suppressWarnings(Mort1Dsmooth(x=x, y=y2, offset=log(e),
                                        w=ww2*wei[w2==1],
                                        method=3, lambda=10^2,
                                        ndx = deg, deg = 3, pord = 2))
  y2.st0 <- exp(XX[[2]] %*% fit2$coef)*e
  y2.st <- rep(0,n)
  y2.st[which(w2==1)] <- y2.st0
  
  ## fitting P-splines for the accident component
  # y1a <- c(y-y2.st) # deleted this line because sometimes y2.st > y => y1a < 0 => error
  y1a <- y[which(w1 == 1)]
  y1a[which(y1a<=0)] <- 0
  y1 <- y1a[which(w1==1)]
  ww1 <- rep(0, length(x))
  ww1[x>3 & x<x1] <- 1
  wei1 <- wei[which(w1==1)]
  fit1 <- suppressWarnings(Mort1Dsmooth(x=x, y=y1, offset=log(e),
                                        w=ww1*wei1,
                                        method=3, lambda=10^4,
                                        ndx = deg, deg = 3, pord = 3))
  y1.st0 <- exp(XX[[1]] %*% fit1$coef)*e
  y1.st <- rep(0,n)
  y1.st[which(w1==1)] <- y1.st0
  
  lhaz1.st <- log(y1.st/e)
  lhaz2.st <- log(y2.st/e)
  y.st <- y1.st + y2.st
  lhaz.st <- log(y.st/e)
  
  ## ## plotting starting values
  ## plot(x, lhaz.act, ylim=c(-12,2))
  ## lines(x, lhaz.st, col=2, lwd=3)
  ## lines(x, lhaz1.st, col=3, lwd=2, lty=3)
  ## lines(x, lhaz2.st, col=4, lwd=2, lty=3)
  
  ## concatenating starting coefficients
  coef.st <- as.vector(c(fit1$coef,
                         fit2$coef))
  
  ## shape penalties
  kappa <- 10^6
  ## including log-concaveness for accident-hump
  D1con <- diff(diag(nb), diff=2)
  w1con <- rep(0, nrow(D1con))
  W1con <- diag(w1con)
  ## including monotonicy for aging part
  D2mon <- diff(diag(nb), diff=1)
  w2mon <- rep(0, nrow(D2mon))
  W2mon <- diag(w2mon)
  
  ## component specific penalty terms
  P1 <- lambda1 * tDD1
  P2 <- lambda2 * tDD2
  ## overall penalty term
  P <- matrix(0,np,np)
  P[1:nb,1:nb] <- P1
  P[1:nb+nb,1:nb+nb] <- P2
  
  coef <- coef.st
  
  ## plotting actual log-mortality
  if(plotIT) plot(x, lhaz.act, ylim=c(floor(min(log(y/e)[y != 0])),0))
  
  Pr <- ridge*diag(np)
  
  ## iterations for a given lambdas-combination
  for(it in 1:max.it){
    
    ds <- 1 - it/max.it
    
    ## penalty for the shape constraints
    P1con <- kappa * t(D1con) %*%W1con%*% D1con
    P2mon <- kappa * t(D2mon) %*%W2mon%*% D2mon
    Psha <- matrix(0,np,np)
    Psha[1:nb,1:nb] <- P1con
    Psha[1:nb+nb,1:nb+nb] <- P2mon
    
    
    ## linear predictor
    eta <- numeric(nx*n)
    for(i in 1:nx){
      eta0 <- rep(0, n)
      eta0[which(indF[,i]==1)] <- XX[[i]] %*% coef[ind[[i]]]
      eta[1:n+(i-1)*n] <- eta0
      if(plotIT){ ## plotting each component
        lines(x[which(indF[,i]==1)],
              eta0[which(indF[,i]==1)], col=i+2)
      }
    }
    ## components
    gamma <- exp(eta)*c(indF)
    ## expected values
    mu <- numeric(n)
    for(i in 1:nx){
      mu <- (e * gamma[1:n+(i-1)*n]) + mu
    }
    ## plotting overall log-mortality
    if(plotIT) lines(x, log(mu/e), col=2, lwd=4)
    ## weights for the IWLS
    w <- mu
    ## modified model matrix for a CLM
    U <- matrix(NA, n, sum(nc))
    for(i in 1:nx){
      u <- gamma[1:n+(i-1)*n]/mu * e
      XXi <- matrix(0, nrow=n, ncol=nc[i])
      XXi[which(indF[,i]==1), ] <- XX[[i]]
      U0 <- u * XXi
      U[,ind[[i]]] <- U0
    }
    U[is.nan(U)] <- 0
    ## regression parts for the P-CLM
    tUWU <- t(U) %*% (w*wei * U)
    tUWUpP <- tUWU + P + Psha + Pr
    r <- y - mu
    tUr <- t(U) %*% r
    ## updating coefficients with a d-step
    coef.old <- coef
    coef <- solve(tUWUpP, tUr + tUWU %*% coef)
    coef <- ds*coef.old + (1-ds)*coef
    ## update weights for shape constraints
    ## accident-hump, log-concaveness
    W1con.old <- W1con
    coef1 <- coef[ind[[1]]]
    W1con <- diag(diff(coef1, diff=2) >= 0)
    ## aging, monotonicity
    W2mon.old <- W2mon
    coef2 <- coef[ind[[2]]]
    W2mon <- diag(diff(coef2) <= 0)
    ## convergence criterion for coefficients
    dif.coef <- max(abs(coef.old-coef))/max(abs(coef))
    ## stopping loop at convergence
    if(dif.coef < 1e-04 & it > 4) break
    ## check convergence
    if(mon) cat(it, dif.coef, "\n")
    ## locator(1)
  }
  ## compute deviance
  yy <- y
  yy[y==0] <- 10^-8
  mumu <- mu
  mumu[mu==0] <- 10^-8
  dev <- 2*sum(y * log(yy/mumu))
  ## effective dimensions
  H <- solve(tUWUpP, tUWU)
  diagH <- diag(H)
  ed1 <- sum(diagH[ind[[1]]])
  ed2 <- sum(diagH[ind[[2]]])
  ## BIC
  bic <- dev + log(n)*sum(diagH)
  
  ## fitted coefficients
  coef.hat <- coef
  ## linear predictor for each component
  etas <- NULL
  for(i in 1:nx){
    etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
  }
  eta1.hat <- etas[[1]]
  eta2.hat <- etas[[2]]
  
  ## linear predictor in a matrix over the whole x
  ETA.hat <- matrix(NA, n, nx)
  for(i in 1:nx){
    ETA.hat[which(indF[,i]==1),i] <- etas[[i]]
  }
  ## linear predictor for overall mortality
  eta.hat <- log(apply(exp(ETA.hat), 1, sum, na.rm=TRUE))
  
  ## fitted values for each component
  gamma1.hat <- exp(eta1.hat)
  gamma2.hat <- exp(eta2.hat)
  
  ## expected values
  mu.hat <- exp(eta.hat)*e
  
  ## deviance residuals
  res1 <- sign(y - mu.hat)
  res2 <- sqrt(2 * (y * log(y/mu.hat) - y + mu.hat))
  res <- res1 * res2
  
  if(plotFIT){
    ## plotting actual and fitted log-mortality
    testofit <- c(expression(paste("ln(",y, "/e)")),
                  expression(paste("ln(", mu, "/e)")),
                  expression(paste("ln(",gamma[1],")")),
                  expression(paste("ln(",gamma[2],")")))
    plot(x, lhaz.act, col=8, pch=16, ylim=c(floor(min(log(y/e)[y!=0])), 0),
         xlab="age",
         ylab="log-mortality",las=1)
    lines(x, eta.hat, col=2, lwd=2)
    lines(x, eta1.hat, col=3)
    lines(x, eta2.hat, col=4)
    legend("top", legend=testofit, col=c(8,2:4),
           lty=c(-1,1,1,1), pch=c(16,-1,-1,-1),ncol=2)
    title(lab)
  }
  ## return object
  out <- list(data=data, coef=coef.hat, XX=XX,
              eta1=eta1.hat, eta2=eta2.hat, eta=eta.hat,
              gamma1=gamma1.hat, gamma2=gamma2.hat,
              mu=mu.hat, resid=res,
              dev=dev, ed1=ed1, ed2=ed2, bic=bic, it=it, maxit=max.it,
              B=B,deg=deg)
  
  return(out)
}

MortSmooth_bbase <-  function(x, xl, xr, ndx, deg){
    ## Input:
    ## x   = abcissae of data
    ## xl  = left boundary
    ## xr  = right boundary
    ## ndx = number of internal knots -1
    ##       or number of internal intervals
    ## deg = degree of the splines
    
    ## Output:
    ## B = matrix with the B-spline basis
    
    ## distance between knots
    dx <- (xr - xl) / ndx
    ## One needs (ndx+1) internal knots and 
    ## deg knots on both right and left side
    ## in order to joint all the B-splines
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    ## Truncated deg-th power functions
    ## equally-spaced at given knots along x
    P <- outer(x, knots, MortSmooth_tpower, deg)
    ## number of B-splines which equal to the number of knots
    n <- dim(P)[2]
    ## in the numerator we have the matrix
    ## of deg+1 differences for each knots
    ## this matrix is rescaled by 
    ## (deg! * dx^deg) == (gamma(deg + 1) * dx ^ deg)
    D <- diff(diag(n), diff = deg + 1) /
      (gamma(deg + 1) * dx ^ deg)
    ## the matrix of differences is used to compute B-splines
    ## as differences of truncated power functions
    ## in P %*% t(D)
    ## the last factor is (-1) ^ (deg + 1)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B
  }

MortSmooth_tpower <-  function(x, t, p){
    ## Input:
    ## x = abcissae of data
    (x - t) ^ p * (x > t)
    ## (x-t)^p gives the curve
    ## (x>t) is an indicator function; it is 1 when x>t
    ## and 0 when x<=t, i.e. before each knot
  }

Mort1Dsmooth <-  function(x, y, offset, w,
           overdispersion=FALSE,
           ndx=floor(length(x)/5),
           deg=3, pord=2, 
           lambda=NULL, df=NULL, method=1,
           coefstart=NULL,
           control=list()){
    ## Input:
    ## x: abcissae of data
    ## y: count response
    ## offset: an a priori known component (optional)
    ## w: a vector of weights to be used in the
    ##    fitting process (optional)
    ## overdispersion: logical on the presence of an
    ##                 overdispersion parameter.
    ##                 Default: FALSE
    ## ndx: number of internal knots -1.
    ##      Default: floor(length(x)/5)
    ## deg: degree of the B-splines. Default: 3
    ## pord: order of differences. Default: 2
    ## lambda: smoothing parameter (optional)
    ## df: a number which specifies the degrees
    ##     of freedom (optional)
    ## method: the method for controlling
    ##         the amount of smoothing. Default: 1
    ## coefstart: eventual starting coefficients
    ## control: a list of control parameters
    ## MON: logical on monitoring
    ## TOL1: convergence of the IWLS algorithm.
    ##       Default: 1e-06.
    ## TOL2: difference between two adjacent
    ##       smoothing parameters in the grid search,
    ##       log-scale. Default: 0.1.
    ## RANGE: range in between lambda is searched, log-scale.
    ##        Default: [10^-8 ; 10^8]
    ## MAX.IT: the maximum number of iterations. Default: 50
    
    
    
    ## Output: a Mort1Dsmooth object containing
    ## aic: Akaike Information Criterion
    ## bic: Bayesian Information Criterion
    ## dev: Poisson-deviance
    ## lev: diag of the hat-matrix
    ## df: effective dimension
    ## lambda: the selected (given) lambda
    ## psi2: (estimated) overdispersion parameter
    ## ndx: number of internal knots -1
    ## deg: degree of the B-splines
    ## pord: order of differences
    ## B: B-splines basis
    ## x: abcissae of data
    ## y: count responses
    ## offset: an a priori known component
    ## w: (only for weighted fits) the specified weights
    ## coef: fitted (penalized) B-splines coefficients  
    ## residuals: the deviance residuals
    ## fitted.values: fitted counts
    ## linear.predictor: fitted linear predictor
    ## logmortality: fitted mortality rates in log-scale
    ## call: the matched call
    ## n: number of observations
    ## tolerance: convergence tolerance
    
    ## checkings:
    if(missing(w)){
      w <- rep(1, length(y))
    }
    ## about x and y:
    if(missing(y)){
      if(is.list(x)){
        if(any(is.na(match(c("x", "y"), names(x)))))
          stop("cannot find x and y in list")
        x <- x$x
        y <- x$y
      }
      else if(is.complex(x)) {
        y <- Im(x)
        x <- Re(x)
      }
      else if(is.matrix(x) && ncol(x) == 2) {
        y <- x[, 2]
        x <- x[, 1]
      }
      else {
        y <- x
        x <- time(x)
      }
    }
    ## about the offset
    m <- length(y)
    if(missing(offset)) {
      offset <- rep(0, m)
    }else{
      if(length(offset) == 1) {
        offset <- rep(offset, m)
      }
      else{
        offset=offset
      }
    }  
    check <- Mort1Dsmooth_checker(x=x, y=y,
                                  offset=offset, w=w,
                                  overdispersion=overdispersion,
                                  ndx=ndx, deg=deg,
                                  pord=pord, 
                                  lambda=lambda, df=df,
                                  method=method,
                                  coefstart=coefstart,
                                  control=control)
    x <- check$x
    y <- check$y
    m <- check$m
    offset <- check$offset
    offsetINIT <- check$offsetINIT
    wei <- check$w
    over <- check$overdispersion
    ndx <- check$ndx
    deg <- check$deg
    pord <- check$pord
    lambda <- check$lambda
    df <- check$df
    MET <- check$method
    a.init <- check$coefstart
    MON <- check$control$MON
    TOL1 <- check$control$TOL1
    TOL2 <- check$control$TOL2
    RANGE <- check$control$RANGE
    MAX.IT <- check$control$MAX.IT
    call <- match.call()
    ## B-splines basis
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
    ## penalty stuff
    nb <- ncol(B)
    D. <- diff(diag(nb), diff=pord)
    DtD <- t(D.)%*%D.
    ## General initialize:
    if(is.null(a.init)){
      y[is.na(y)] <- 0
      ## 0) simple poisson-GLM with ages(years)
      ##    only for the interpolation cases
      fit0 <- glm(round(y) ~ x + offset(offset),
                  family=poisson, weights=wei)
      ## 1) simple penalized-poisson-GLM with B
      etaGLM <- log(fit0$fitted) - offset
      eta0 <- log((y + 1)) - offset
      eta0[wei==0] <- etaGLM[wei==0]
      mu0 <- exp(offset + eta0)
      w0 <- wei*mu0
      z0 <- wei*((y - mu0)/mu0 + eta0)
      BtWB <- t(B) %*% (w0 * B)
      BtWz <- t(B) %*% (w0 * z0)
      a.init <- solve(BtWB + 1e08 * DtD, BtWz)
    }else{
      a.init=a.init
    }
    psi2 <- 1
    
    ## ## plotting starting coef
    ## ran <- range(log(y/e), B%*%a.init,
    ##              na.rm=TRUE, finite = TRUE)
    ## plot(x, log(y/e), ylim=ran)
    ## lines(x, B%*%a.init, col=2, lwd=2)
    
    ## optimize AIC or BIC
    if(MET==1|MET==2){
      ## if overdisperion is true
      if(over){
        tol.over <- 10
        i.over <- 0
        while(tol.over > 1e-03 && i.over < 5){
          i.over <- i.over+1
          lambda.hat <- Mort1Dsmooth_optimize(x=x,
                                              y=y,
                                              offset=offset,
                                              wei=wei,
                                              psi2=psi2,
                                              B=B, DtD=DtD,
                                              a.init=a.init,
                                              MON=MON,
                                              TOL1=TOL1,
                                              TOL2=TOL2,
                                              RANGE=RANGE,
                                              MAX.IT=MAX.IT,
                                              MET=MET)
          FIT <- Mort1Dsmooth_estimate(x=x, y=y,
                                       offset=offset,
                                       wei=wei,
                                       psi2=psi2,
                                       B=B, 
                                       lambda=lambda.hat,
                                       DtD=DtD,
                                       a.init=a.init, 
                                       MON=MON, TOL1=TOL1,
                                       MAX.IT=MAX.IT)
          ## recalculating overdispersion parameter
          psi2.old <- psi2
          psi2 <- FIT$dev / (m - FIT$df)
          tol.over <- abs(psi2 - psi2.old)/abs(psi2)
        }
      }else{## if psi2==1
        lambda.hat <- Mort1Dsmooth_optimize(x=x,
                                            y=y,
                                            offset=offset,
                                            wei=wei,
                                            psi2=psi2,
                                            B=B, DtD=DtD,
                                            a.init=a.init,
                                            MON=MON,
                                            TOL1=TOL1,
                                            TOL2=TOL2,
                                            RANGE=RANGE,
                                            MAX.IT=MAX.IT,
                                            MET=MET)
        FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                     wei=wei,
                                     psi2=psi2,
                                     B=B, 
                                     lambda=lambda.hat,
                                     DtD=DtD,
                                     a.init=a.init, 
                                     MON=MON, TOL1=TOL1,
                                     MAX.IT=MAX.IT)
        psi2 <- FIT$dev / (m - FIT$df)
      }
      if(log10(lambda.hat)>=log10(RANGE[2]) |
           log10(lambda.hat)<=log10(RANGE[1])) {
        warning(paste("optimal lambda at the edge of the grid."))
      }
    }
    ## given lambda
    if(MET==3){
      lambda.hat <- lambda
      FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                   wei=wei,
                                   psi2=psi2,
                                   B=B, 
                                   lambda=lambda.hat, DtD=DtD,
                                   a.init=a.init, 
                                   MON=MON, TOL1=TOL1,
                                   MAX.IT=MAX.IT)
      psi2 <- FIT$dev / (m - FIT$df)
    }
    ## optimize given df
    if(MET==4){
      Mort1Dsmooth_opt_df <- function(X){
        FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                     wei=wei,
                                     psi2=psi2,
                                     B=B, 
                                     lambda=X, DtD=DtD,
                                     a.init=a.init, 
                                     MON=MON, TOL1=TOL1,
                                     MAX.IT=MAX.IT)
        return(abs(FIT$df - df))
      }
      by.lambda <- length(seq(log10(RANGE[1]),
                              log10(RANGE[2]),by=TOL2))
      lambda.hat <- cleversearch(fn=Mort1Dsmooth_opt_df,
                                 lower=log10(RANGE[1]),
                                 upper=log10(RANGE[2]),
                                 ngrid=by.lambda,
                                 startvalue=1,
                                 logscale=TRUE,
                                 verbose=FALSE)[[1]]
      if(log10(lambda.hat)>=log10(RANGE[2]) |
           log10(lambda.hat)<=log10(RANGE[1])){
        warning(paste("optimal lambda at the edge of the grid."))
      }
      FIT <- Mort1Dsmooth_estimate(x=x, y=y, offset=offset,
                                   wei=wei,
                                   psi2=psi2,
                                   B=B, 
                                   lambda=lambda.hat,
                                   DtD=DtD,
                                   a.init=a.init, 
                                   MON=MON, TOL1=TOL1,
                                   MAX.IT=MAX.IT)
      psi2 <- FIT$dev / (m - FIT$df)
    }
    aic <- FIT$aic
    bic <- FIT$bic
    df <- FIT$df
    dev <- FIT$dev
    coef <- FIT$a
    psi2 <- psi2
    h <- FIT$h
    tolerance <- FIT$tol
    eta.hat <- B%*%coef
    logmortality <- eta.hat
    fitted.values <- exp(eta.hat + offset)
    res <- sign(y - fitted.values) *
      sqrt(2 * (y * log(ifelse(y == 0, 1,
                               y/fitted.values)) -
                  (y - fitted.values)))
    lin <- as.vector(eta.hat) + as.vector(offset)
    ## output
    object <- list(## fitted values
      coefficients=as.vector(coef),
      residuals=as.vector(res),
      fitted.values=as.vector(fitted.values),
      linear.predictors=lin,
      logmortality=logmortality,
      ## diagnostics
      lev=h, df=df, deviance=dev,
      aic=aic, bic=bic,
      psi2=psi2, lambda=lambda.hat,
      ## general
      call=call, n=m, tolerance=tolerance,
      ## smoothing specifications
      ndx=ndx, deg=deg, pord=pord, B=B,
      ## overdispersion=over,
      ## observed values
      x=x, y=as.vector(y),
      offset=as.vector(offsetINIT),
      w=as.vector(wei)
    )
    class(object) <- "Mort1Dsmooth"
    object
  }

Mort1Dsmooth_checker <-
  function(x, y, offset, w,
           overdispersion,
           ndx, deg, pord, 
           lambda, df, method,
           coefstart,
           control){
    ## Input:
    ## x: abcissae of data
    ## y: count response
    ## offset: an a priori known component (optional)
    ## w: weights
    ## overdispersion: logical on the presence of
    ##                  overdispersion
    ## ndx: number of internal knots -1.
    ##      Default: floor(length(x)/5)
    ## deg: degree of the B-splines. Default: 3
    ## pord: order of differences. Default: 2
    ## lambda: smoothing parameter. Default: NULL (optional)
    ## df: a number which specifies the degrees of freedom.
    ##     Default: NULL (optional)
    ## method: the method for controlling the amount of
    ##         smoothing. Default: 4
    ## coefstart: eventual initial coefficients
    ## control: a list of control parameters
    
    ## Output: a list containing CHECKED arguments
    ##         for the Mort1Dsmooth function
    m <- length(y)
    
    offsetINIT <- offset
    ## about infinitive or NA offset
    whioff <- which(is.infinite(offset) | is.na(offset))
    whiwei <- which(w==0)
    if(any(!whioff%in%whiwei)){
      stop("weights different from zero associated with infinitive or NA offset values")
    }
    offset[c(whioff, whiwei)] <- 100
    ## about lengths and wrong values
    if (length(x)!=length(y)) 
      stop("Arguments must have same length")
    if (length(y) != m | length(offset) != m) 
      stop("Argument arrays of wrong length")
    if (deg < 1 | deg >= 10) 
      stop("Wrong value for deg")
    if (pord <= 0 | pord >= 5) 
      stop("Wrong value for pord")
    if (ndx < 2 | ndx >= floor(m*.9)) 
      stop("Wrong value for ndx")
    coefstart.check <- is.null(coefstart)
    if(!coefstart.check){
      if(length(coefstart)!=(ndx+deg)){
        stop("coefstart must have length equal to ndx+deg")
      }
    }
    ## about method
    if (method != 1 & method != 2 &
          method != 3 & method != 4) 
      stop("Wrong value for method")
    ## method = 1 adjusts lambda so that the BIC is minimized
    ## method = 2 adjusts lambda so that the AIC is minimized
    ## method = 3 uses the value supplied for lambda
    ## method = 4 adjusts lambda so that the degrees of
    ##          freedom is equal to df
    ## check-point methods
    lambda.check <- is.null(lambda)
    df.check <- is.null(df)
    MET <- NULL
    ## both lambda and df NULL
    if(lambda.check & df.check & method==1){MET=1}
    if(lambda.check & df.check & method==2){MET=2}
    if(lambda.check & df.check & method==3){
      stop("with method 3, provide lambda")
    }
    if(lambda.check & df.check & method==4){
      stop("with method 4, provide df")
    }
    ## lambda NULL and df GIVEN
    if(lambda.check & !df.check & method==1){
      stop("df and method 1 cannot be chosen together")
    }
    if(lambda.check & !df.check & method==2){
      stop("df and method 2 cannot be chosen together")
    }
    if(lambda.check & !df.check & method==3){
      stop("df and method 3 cannot be chosen together")
    }
    if(lambda.check & !df.check & method==4){MET=4}
    ## lambda GIVEN and df NULL
    if(!lambda.check & df.check & method==1){
      stop("lambda and method 1 cannot be chosen together")
    }
    if(!lambda.check & df.check & method==2){
      stop("lambda and method 2 cannot be chosen together")
    }
    if(!lambda.check & df.check & method==3){MET=3}
    if(!lambda.check & df.check & method==4){
      stop("lambda and method 4 cannot be chosen together")
    }
    ## both lambda and df GIVEN, regardless method
    if(!lambda.check & !df.check){
      stop("lambda and df cannot be chosen together")
    }  
    ## impossible values for lambda and df
    if(!lambda.check && lambda<0)
      stop("lambda must be positive")
    if(!df.check && df<pord)
      stop("df must be larger than pord")
    if(!df.check && df>c(ndx+deg))
      stop("df must be smaller than ndx+deg")
    if (!df.check & length(df)!=1)
      stop("df must be length 1")
    if (!lambda.check & length(lambda)!=1)
      stop("lambda must be length 1")
    ## setting control-parameters
    con <- list(MON=FALSE, TOL1=1e-06, TOL2=0.5,
                RANGE=c(10^-4, 10^6), MAX.IT=50)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
      warning("unknown names in control: ",
              paste(noNms, collapse = ", "))
    ## stop about weights
    if(length(w)!=m){ 
      stop("length of w and y must be equal")
    }
    ## warning about interpolation/extrapolation
    if(any(w==0)){
      warning("Interpolation and/or extrapolation is taking place", call. = FALSE)
    }
    ## about the overdispersion parameter
    if(overdispersion & method==3)
      warning("given method 3, overdispersion is computed a posteriori")
    if(overdispersion & method==4)
      warning("given method 4, overdispersion is computed a posteriori")
    ## warning about weights
    if(min(w) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    ## returning
    llist <- list(x=x, y=y, offset=offset, w=w,
                  offsetINIT=offsetINIT,
                  overdispersion=overdispersion, m=m,
                  ndx=ndx, deg=deg, pord=pord, 
                  lambda=lambda, df=df, method=method,
                  coefstart=coefstart,
                  control=con)
    llist
  }

HMD2MH <- function(country, year, dim = "period", xtra = FALSE, sex, min = 0, max = NULL, username = NULL, password = NULL, path = NULL){
  
  if(all(is.null(c(username,password,path)))){warning("You must provide either valid username and password (web) or path to the folder (local).")}
  
  if(sex %in% c("males","females","total") == F){warning("sex must be one of 'females', 'males' or 'total'")}
  
  if(dim %in% c("period","cohort") == F){warning("dim must be one of 'period' or 'cohort'")}
  
  if(dim == "period"){
    
    if(is.null(path) == F){
      
      dx <- readHMD(filepath = file.path(path,country,"STATS","Deaths_1x1.txt"), fixup = T)
      nx <- readHMD(filepath = file.path(path,country,"STATS","Exposures_1x1.txt"), fixup = T)
      
    }
    
    if(all(is.null(c(username,password)) == F)){
      
      dx <- readHMDweb(username = username, password = password, CNTRY = country, item = "Deaths_1x1", fixup = T)
      nx <- readHMDweb(username = username, password = password, CNTRY = country, item = "Exposures_1x1", fixup = T)
      
    }
    
    dx <- dx[dx$Year == year,]
    nx <- nx[nx$Year == year,]
    
    if(is.null(max)){max <- max(dx$Age)}
    
    dx <- dx[dx$Age <= max,]
    nx <- nx[nx$Age <= max,]
    
    x <- dx$Age
    
    if(sex == "females"){dx <- dx[,3] ; nx <- nx[,3]}
    if(sex == "males"){dx <- dx[,4]; nx <- nx[,4]}
    if(sex == "total"){dx <- dx[,5]; nx <- nx[,5]}
    
    data <- data.frame(x = x, d = dx, n = nx)
    
    data$m <- dx / nx
    
  }
  
  if(dim == "cohort"){
    
    if(xtra == TRUE & is.null(path)){warning("The extrapolation option is only compatible with locally stored HMD data.")}
    
    cnx <- readHMD(filepath = paste(path,country,"STATS","cExposures_1x1.txt",sep = "/"), fixup = T)
    cmx <- readHMD(filepath = paste(path,country,"STATS","cMx_1x1.txt",sep = "/"), fixup = T)
    
    cnx <- cnx[cnx$Age <= max,]
    cmx <- cmx[cmx$Age <= max,]
    
    if(xtra == TRUE){cmx <- xtp(path = path, max = max, country = country, year = year, sex = sex, cmx = cmx)}else{
      cmx <- cmx[cmx$Year == year,]
    }
    cnx <- cnx[cnx$Year == year,]
    
    x <- cmx$Age
    
    if(sex == "females"){cmx <- cmx[,3] ; cnx <- cnx[,3]}
    if(sex == "males"){cmx <- cmx[,4]; cnx <- cnx[,4]}
    if(sex == "total"){cmx <- cmx[,5]; cnx <- cnx[,5]}
    
    if(xtra == TRUE){ for(i in which(is.na(cnx))){cnx[i] <- cnx[i-1] - cnx[i-1] * cmx[i]}  }
    
    cdx <- round(cmx * cnx)
    cmx <- cdx / cnx
    
    data <- data.frame(x = x, d = cdx, n = cnx, m = cmx)
    
  }
  
  return(data)
  
}

HCD2MH <- function(country, year, sex, list, unabr = FALSE, path){
  
  if(sex %in% c("males","females","total") == F){warning("sex must be one of 'females', 'males' or 'total'")}
  
  mxc <- read.csv(file = file.path(path,country,paste(country,"_m_",list,"_idr.csv",sep = "")))
  dxc <- read.csv(file = file.path(path,country,paste(country,"_d_",list,"_idr.csv",sep = "")))
  nx  <- read.csv(file = file.path(path,country,paste(country,"_e.csv",sep = "")))
  
  if(sex == "females"){sex <- 2}
  if(sex == "males"){sex <- 1}
  if(sex == "total"){sex <- 3}
  
  mxc <- mxc[mxc$year == year & mxc$sex == sex,]
  dxc <- dxc[dxc$year == year & dxc$sex == sex,]
  nx  <- nx[nx$year == year & nx$sex == sex,]
  
  af <- unique(c(mxc$agf,dxc$agf))
  if(length(af) > 1){warning("Two different age formats are present in the selected dataset.")}
  #data("agf")
  age <- as.character(agf[,af])
  
  if(0 %in% mxc$cause){mxc <- mxc[mxc$cause != 0,] ; dxc <- dxc[dxc$cause != 0,]}
  causes <- mxc$cause
  
  mxc <- as.data.frame(t(mxc[,7:ncol(mxc)]))
  suppressWarnings(mxc <- as.data.frame(apply(mxc, 2, as.numeric)))
  mxc <- mxc / 1e6
  mxc <- mxc[rowSums(is.na(mxc)) < ncol(mxc),]
  if(af == 2){mxc <- mxc[-19,]}
  if(af == 3){mxc <- mxc[-c(19,21),]}
  if(af == 4){mxc <- mxc[-c(19,21,23),]}
  names(mxc) <- causes
  mxc <- mxc[,colSums(is.na(mxc)) < nrow(mxc)]
  
  
  dxc <- as.data.frame(t(dxc[,8:ncol(dxc)]))
  suppressWarnings(dxc <- as.data.frame(apply(dxc, 2, as.numeric)))
  dxc <- dxc[rowSums(is.na(dxc)) < ncol(dxc),]
  if(af == 2){dxc <- dxc[-19,]}
  if(af == 3){dxc <- dxc[-c(19,21),]}
  if(af == 4){dxc <- dxc[-c(19,21,23),]}
  names(dxc) <- causes
  dxc <- dxc[,colSums(is.na(dxc)) < nrow(dxc)]
  
  
  nx <- t(nx[1,6:30])
  if(af == 1){nx <- nx[-c(20:length(nx))]}
  if(af == 2){nx <- nx[-c(19,22:length(nx))]}
  if(af == 3){nx <- nx[-c(19,21,24:length(nx))]}
  if(af == 4){nx <- nx[-c(19,21,23)]}
  
  if(af == 1){x <- age[-c(20:length(age))]}
  if(af == 2){x <- age[-c(19,22:length(age))]}
  if(af == 3){x <- age[-c(19,21,24:length(age))]}
  if(af == 4){x <- age[-c(19,21,23)]}
  age <- x
  plus <- grep(x = as.character(x), pattern = "+", fixed = TRUE)
  x[plus] <- sub(x = x[plus], pattern = "+", replacement = "", fixed = TRUE)
  inter <- do.call(rbind,strsplit(as.character(x), "-"))
  inter <- cbind(apply(inter,2,as.numeric))
  inter[,2] <- inter[,2] + 1
  inter[nrow(inter),2] <- inter[nrow(inter),1] + 10
  x <- rowMeans(inter)
  
  if(unabr == TRUE){
    
    ua <- unabridge(dxc = dxc, nx = nx, inter = inter)
    age <- ua$x
    x <- ua$x + 0.5
    dxc <- ua$dxc
    nx <- ua$nx
    mxc <- dxc / nx
    inter <- cbind(x[-length(x)],x[-1])
    
  }
  
  #data("lists")
  if(list == "short"){lab <- as.data.frame(short)}
  if(list == "interm"){lab <- as.data.frame(interm)}
  if(list == "full"){message("The labels of the full list must be downloaded from the HCD website as they differ by country.") ; lab <- NA}
  lab <- lab[match(names(mxc),lab$no),]
  
  return(list(mxc = mxc, dxc = dxc, nx = nx, x = x, age = age, inter = inter, lab = lab))
  
}

.onAttach <- function(libname, pkgname){
  suppressWarnings(descr <- utils::packageDescription("MortHump"))
  # 	if(utils::packageVersion("MortHump")$minor %% 2 == 0) {
  # 		state <- "stable"
  # 	}
  # 	else {
  #  		state <- "development"
  #  	}
  state <- ""
  if(!is.null(descr$Built)){
    builtDate <- paste(" (Built: ", strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2], ")", sep="")
  }else{
    builtDate <- ""
  }
  packageStartupMessage("This is MortHump ", state, " version ", descr$Version, builtDate)
  packageStartupMessage("\nTo cite MortHump in publications please use:")
  packageStartupMessage("Remund, Camarda and Riffe. Analyzing the Young Adult Mortality Hump in R with MortHump.")
  packageStartupMessage("MPIDR Technical Report TR-2018-3. (2018).")
}

summary.morthump <- function(object,...){
  
  fit <- object
  
  model <- attr(fit, "model")
  data <-  fit$data
  x <- data$x
  if(x[1] == 0){
    if(model %in% c("hp","hpk","hps")){x[1] <- 1e-5}else{x[1] <- 0.9}
  }
  
  
  # mx fitted, hump, and w/o hump
  mhat <- predict(object = fit, x = x, which = "fitted")
  mhump <- predict(object = fit, x = x, which = "hump")
  mhyp <- mhat - mhump
  
  par1 <- coef(object = fit)
  
  if(model %in% c("hp","hpk")){
    
    # Fisher test compared to a model w/o hump
    w <- fit$w
    m   <- data$m
    fit2 <- hps.fit(data = data, w = w)
    par2 <- coef(object = fit2)
    mhat2 <- predict(object = fit2)
    df1 <- length(par1) - length(par2)
    df2 <- length(x) - length(par)
    res1 <- sum(w * (m - mhat)^2, na.rm = TRUE)
    res2 <- sum(w * (m - mhat2)^2, na.rm = TRUE)
    test <- ((res2 - res1) / df1) / (res2 / df2)
    suppressWarnings(pval <- pf(test, df1 = df1, df2 = df2, lower.tail = FALSE))
    #pval <- stats:::anovalist.nls(fit,fit2)$"Pr(>F)"[2]
    
  }else{pval <- NA}
  
  # loss of life expectancy
  e0hat <- LT(Mx = mhat, ages = data$x, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hyp <- LT(Mx = mhyp, ages = data$x, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hump  <- round(e0hyp - e0hat,2)
  
  # loss of lives
  dh <- round(sum(data$n * mhump))
  
  # years of life lost
  yll <- round((data$n * mhump) %*% LT(Mx = mhat, ages = data$x, mxsmooth = FALSE, verbose = FALSE)$ex)
  
  if(model != "hps"){
    
    # functions
    pdf <- function(t){predict(object = fit, x = t, which = "hump", tol = 2) / sum(mhump)}
    pdf <- Vectorize(pdf)
    cdf <- Vectorize(FUN = function(t){integrate(f = function(a){pdf(a)}, lower = min(x), upper = t)$value})
    if(model == "sse"){cdf <- function(t){sum(pdf(seq(min(x),t,0.1))) / 10}}
    cdf <- Vectorize(cdf)
    qtl <- function(p){optimize(f = function(t){abs(cdf(t) - p)}, interval = range(x))$minimum}
    qtl <- Vectorize(qtl)
    
    # descriptive statistics
    mode.st <- x[which.max(predict(object = fit, which = "hump"))]
    mode <- optimize(f = function(x){-pdf(x)}, lower = mode.st * 0.8, upper = mode.st * 1.2)$minimum
    median <- qtl(0.5)
    mean <- as.numeric(pdf(x) %*% x)
    sd <- as.numeric(sqrt(pdf(x) %*% (data$x - mean)^2))
  }
  
  
  
  if(model != "hps"){
    out <- list(mhat = mhat, mhump = mhump, pval = pval, pdf = pdf, cdf = cdf, qtl = qtl, loss = e0hump, dh = dh, yll = yll, mode = mode, median = median, mean = mean, sd = sd, par = par1)
  }else{
    out <- list(mhat = mhat, mhump = mhump, pval = pval, loss = e0hump, dh = dh, yll = yll, par = par1, it = fit$it)
  }
  
  class(out) <- "summary.morthump"
  attr(out, "extra") <- list(model = model, type = fit$type, method = fit$method, it = fit$it, maxit = fit$maxit)
  
  invisible(out)
}

print.summary.morthump <- function(x, ...){
  
  model <- attr(x, "extra")$model
  type <- attr(x, "extra")$type
  
  cat("[[>>]] Information about the model used to fit the mortality schedule")
  
  
  cat("\nType of the model:",type)
  cat("\nName of the model:",model)
  cat("\nAlgorithm used for optimization:",attr(x,"extra")$method)
  if(type == "parametric"){
    cat("\nEstimated coeffients: ", paste0(names(x$par),"=",signif(x$par,3)))
    cat("\nF-test against a model w/o hump: ", signif(x$pval,3),"(H0: this model does not fit better than a model without a hump)")
  }else{
    cat("\nNumber of iterations:",attr(x, "extra")$it,"( max =",attr(x, "extra")$maxit,")")
  }
  
  cat("\n\n")
  
  if(model != "hps"){
    
    cat("[[>>]] Descriptive statistics about the hump")
    
    cat("\nIntensity")
    cat("\nLife expectancy at birth lost due to the hump: ", x$loss, "year")
    cat("\nDeaths due to the hump: ", x$dh)
    cat("\nYears of life lost due to the hump: ", x$yll)
    cat("\n")
    
    cat("\nCentrality")
    cat("\nMode:", round(x$mode,2))
    cat("\nMedian:", round(x$median,2))
    cat("\nMean:", round(x$mean,2))
    cat("\n")
    
    cat("\nDispersion")
    cat("\nStandard deviation:", round(x$sd,2))
    
  }
  
  invisible(x)
}

hp <- function(x, par){
  
  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  E <- par[5]
  F <- par[6]
  G <- par[7]
  H <- par[8]
  
  y <- A ^ ((x + B) ^ C) + D * exp(-E * (log(x) - log(F)) ^ 2) + (G * (H ^ x))/(1 + G * (H ^ x))
  
  names(y) <- NULL
  
  return(y)
}


hp.fit <- function(data, method = "port", w = 1/data$m, start = NULL){
  
  if(data$x[1] == 0){data$x[1] <- 1e-5}
  
  form <- as.formula(m ~ A ^ (x + B) ^ C + D * exp(-E * (log(x) - log(F)) ^ 2) + G * H ^ x / (1 + G * (H ^ x)))
  
  if(is.null(start)){
    
    start <- list(A = 0.001, B = 0.005, C = 0.11, D = 0.0015, E = 8, F = 20, G = 0.00003, H = 1.105)
    
    lower <- c(0.0001, 0.000001, 0.0001, 0, 1, 16, 0.0000001, 0.5)
    
    upper <- c(0.1, 0.5, 1, 0.01, 50, 25, 0.01, 1.5)
    
  }else{
    lower <- start$lower
    upper <- start$upper
    start <- start$start
  }
  
  if(method %in% c("port","lm","gnm","bayes") == FALSE){warning("The method must be one of 'port, 'lm' or 'gnm'.")}
  
  if(method == "port"){
    
    fit <- nls(formula = form, data = data, start = start, lower = lower, upper = upper,
               algorithm = "port", weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "lm"){
    
    fit <- nlsLM(formula = form, data = data, start = start, lower = lower, upper = upper,
                 weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "gnm"){
    
    # see Currie 2014 and Debon et al. 2005
    
    warning("Sorry, the estimation of the Heligman-Pollard model as a generalized linear model is not implemented yet.")
    
  }
  
  if(method == "bayes"){
    
    # fun <- function(m){rnorm(n = 8e3, mean = m, sd = m/10)}
    
    # prior <- do.call(cbind,lapply(start,fun))
    
    # fit <- hp.bm.imis(prior = prior, nrisk = data$exp, ndeath = data$d, K = 10)
    
    warning("Sorry, the estimation of the Heligman-Pollard model using Bayesian statistics is not implemented yet.")
    
  }
  
  if(data$x[1] < 1){data$x[1] <- 0}
  
  fit$data <- data
  fit$method <- method
  fit$w <- w
  fit$type <- "parametric"
  
  return(fit)
  
}

hpk <- function(x, par){
  
  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  E <- par[5]
  F <- par[6]
  G <- par[7]
  H <- par[8]
  k <- par[9]
  
  y <- A ^ (x + B) ^ C + D * exp(-E * (x <= F) * (log(x) - log(F)) ^ 2) * exp(-k * E * (x > F) * (log(x) - log(F)) ^ 2) + G * H ^ x / (1 + G * (H ^ x))
  
  return(y)
}

hpk.fit <- function(data, method = "port", w = 1/data$m, start = NULL){
  
  # according to Heligman and Pollard (1980), weights should be 1/q
  # according to Brillinger (1986), weights should be 1/d
  # according to Heligman and Pollard (1980), the response variable should be q or q/(1-q), but it is defined here as m
  
  
  if(data$x[1] == 0){data$x[1] <- 1e-5}
  
  form <- as.formula(m ~ A ^ (x + B) ^ C + D * exp(-E * (x <= F) * (log(x) - log(F)) ^ 2) * exp(-k * E * (x > F) * (log(x) - log(F)) ^ 2) + G * H ^ x / (1 + G * (H ^ x)))
  
  if(is.null(start)){
    
    start <- list(A = 0.001, B = 0.005, C = 0.11, D = 0.0015, E = 8, F = 20, G = 0.00003, H = 1.105, k = 0.5)
    
    lower <- c(0.0001, 0.000001, 0.0001, 0, 1, 16, 0.0000001, 0.5, 0)
    
    upper <- c(0.1, 0.5, 1, 0.01, 50, 22, 0.01, 1.5, 1)
    
  }else{
    lower <- start$lower
    upper <- start$upper
    start <- start$start
  }
  
  if(method %in% c("port","lm","gnm","bayes") == F){warning("The method must be one of 'port, 'lm' or 'gnm'.")}
  
  if(method == "port"){
    
    fit <- nls(formula = form, data = data, start = start, lower = lower, upper = upper,
               algorithm = "port", weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "lm"){
    
    fit <- nlsLM(formula = form, data = data, start = start, lower = lower, upper = upper,
                 weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "gnm"){
    
    # see Currie 2014 and Debon et al. 2005
    
    warning("Sorry, the estimation of the Kostaki model as a generalized linear model is not implemented yet.")
    
  }
  
  if(method == "bayes"){
    
    # fun <- function(m){rnorm(n = 8e3, mean = m, sd = m/10)}
    
    # prior <- do.call(cbind,lapply(start,fun))
    
    # fit <- hp.bm.imis(prior = prior, nrisk = data$exp, ndeath = data$d, K = 10)
    
    warning("Sorry, the estimation of the Kostaki model using Bayesian statistics is not implemented yet.")
    
  }
  
  if(data$x[1] < 1){data$x[1] <- 0}
  
  fit$data <- data
  fit$method <- method
  fit$w <- w
  fit$type <- "parametric"
  
  
  return(fit)
  
}

coef.morthump <- function(object,...){
  
  fit <- object
  
  par <- unlist(fit$coef)
  
  return(par)
  
}

plot.morthump <- function(x, which,...){
  
  fit <- x
  
  model <- attr(fit, "model")
  data <- fit$data
  if(data$x[1] == 0 & model != "sse"){data$x[1] <- 1e-5}
  
  hump <- summary(fit, print = FALSE)
  
  if(which == 1){
    
    plot(x = data$x, y = log(data$m), type = "p", ylab = "mx (log scale)", xlab = "age", pch = 16, col = "grey", axes = FALSE, ...)
    axis(side = 1)
    yl <- 10^(floor(log10(min(data$m))):0)
    axis(side = 2, at = log(yl), labels = yl, las = 1)
    box()
    title("Age-specific death rates")
    
    lines(x = data$x, y = log(hump$mhat), lwd = 2, col = 1)
    
    if(model != "hps"){
      
      lines(x = data$x, y = log(hump$mhump), lwd = 1, lty = 2)
      
      polygon(x = c(data$x, rev(data$x)), y = log(c(hump$mhat - hump$mhump, rev(hump$mhat))), col = "#00000020")
      
      # text(labels = paste(round(hump$loss,2),"yr"), x = which.max(hump$mhump),y = log(predict(fit, x = which.max(hump$mhump)) / 2))
      
      legend(x = "bottomright", legend = c("observed","fitted","hump"), col = c(8,1,1), pch = c(16,NA,NA), lwd = c(NA,2,1), lty = c(NA,1,2))
      
    }else{legend(x = "bottomright", legend = c("observed","fitted"), col = c(8,1), pch = c(16,NA), lwd = c(NA,2), lty = c(NA,1))}
    
  }
  
  if(which == 2){
    
    if(model == "hps"){stop("The HPS model does not assume any hump")}
    
    mode <- hump$mode
    median <- hump$median
    mean <- hump$mean
    pdf <- hump$pdf
    cdf <- hump$cdf
    qtl <- hump$qtl
    cols <- brewer.pal(n = 3, name = "Set1")
    
    plot(x = data$x, y = pdf(data$x), type = "n", las = 1, xlab = "age", ylab = "pdf", xlim = range(which(pdf(data$x) > 1e-5)), ...)
    
    xx <- seq(min(data$x), max(data$x), 0.1)
    polygon(x = c(xx,max(xx),min(xx)),  y = c(pdf(xx),0,0), col = "grey", border = NA)
    
    # deal with overlapping labels
    xydist <- function(x1, x2){sqrt((x1 - x2)^2 + (pdf(x1) - pdf(x2))^2)}
    d1 <- xydist(mode, median)
    d2 <- xydist(median, mean)
    
    points(x = mode, y = pdf(mode), pch = 3, lwd = 2, col = cols[1])
    text(x = mode,  y = pdf(mode), labels = "mode", pos = ifelse(d1 < 1 & mode < median,2,4), col = cols[1])
    segments(x0 = mode, x1 = mode, y0 = 0, y1 = pdf(mode), lty = 3, col = cols[1])
    
    points(median, pdf(median), pch = 3, lwd = 2, col = cols[2])
    text(x = median,  y = ifelse(d2 < 1 & pdf(median) > pdf(mean), pdf(median)*1.02, pdf(median)), labels = "median", pos = ifelse(median > mode,4,2), col = cols[2])
    segments(x0 = median, x1 = median, y0 = 0, y1 = pdf(median), lty = 3, col = cols[2])
    
    points(mean, pdf(mean), pch = 3, lwd = 2, col = cols[3])
    text(x = mean,  y = ifelse(d2 < 1 & pdf(median) < pdf(mean), pdf(mean)*1.02, pdf(mean)), labels = "mean", pos = ifelse(mean > mode,4,2), col = cols[3])
    segments(x0 = mean, x1 = mean, y0 = 0, y1 = pdf(mean), lty = 3, col = cols[3])
    
    # legend("topright", legend = c(expression(paste(e[0]^-H - e[0],"=",round(hump$loss,2),"yr")),expression(paste(d^H,"=", hump$dh))), box.lty = 0)
    
    title("Density of the hump")
    
  }
  
}

confint.morthump <- function(object, parm = "e0", level = 0.95, method, iter = 1000,...){
  
  fit <- object
  loss <- summary(fit)$loss
  
  # function to shift a variable
  shift <- function(x,shift_by){
    #http://www.r-bloggers.com/generating-a-laglead-variables/
    stopifnot(is.numeric(shift_by))
    stopifnot(is.numeric(x))
    
    if (length(shift_by)>1)
      return(sapply(shift_by,shift, x=x))
    
    out<-NULL
    abs_shift_by=abs(shift_by)
    if (shift_by > 0 )
      out <- c(tail(x, -abs_shift_by), rep(NA, abs_shift_by))
    else if (shift_by < 0 )
      out <- c(rep(NA, abs_shift_by), head(x, -abs_shift_by))
    else
      out <- x
    out
  }
  
  # quick function to compute life expectancy (quicker than LT)
  LE <- function(mx){
    if (is.null(dim(mx))){
      mx <- as.matrix(mx)
    }
    N             <- nrow(mx)
    ax            <- mx * 0 + .5
    ax[1, ]       <- 0.1
    qx            <- mx / (1 + (1 - ax) * mx)
    qx[N, ]       <- 1
    qx[qx > 1]    <- 1
    px            <- 1 - qx
    lx            <- apply(px, 2,function(.px,.N){
      c(1,cumprod(.px)[-.N])
    }, .N = N)
    dx            <- lx * qx
    dx[N, ]       <- 1 - colSums(dx[-N, , drop = FALSE])
    Lx            <- rbind(
      lx[2:N, , drop = FALSE] + ax[1:(N -  1), , drop = FALSE] * dx[1:(N - 1), , drop = FALSE],
      mx[N,] / ax[N,]
    )
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] <- 0
    Tx            <- apply(Lx, 2,function(.Lx){
      rev(cumsum(rev(.Lx)))
    })
    ex            <- Tx / lx
    return(ex[1, ])
  }
  
  data <- fit$data
  
  if(parm != "e0"){warning("Confidence intervals are only available for e0.")}
  if(sum(data$n) < 5000){message("Population under exposure lower than 5000. The variance estimates will not be reliable (see Eayres & William 2004 and Toson & Baker 2003:10).")}
  if(any(data$d == 0)){message("At least one age group has zero deaths. This may harm the computation of the variance, but the effect for populations over 5000 is minimal (see Toson & Baker 2003:17).")}
  
  stopifnot(method %in% c("MC", "chiang"))
  
  Mx    <- summary(fit, print = FALSE)$mhat
  n     <- length(Mx)
  Nx    <- round(data$n)
  Dx    <- Nx * Mx
  e0    <- LE(Mx)
  
  if(method == "MC"){
    rmx         <- t(mapply(function(.Mx, .N, iter){
      rbinom(n = iter, prob = .Mx, size = .N)
    },Mx,Nx,iter=iter)) / Nx
    re0         <- LE(rmx)
    lower       <- quantile(re0, probs = (1 - level) / 2)
    upper       <- quantile(re0, probs = 1 - (1 - level) / 2)
  }
  
  if(method == "chiang"){
    if(any(Dx < 5)){
      message("Number of deaths is too low in at least one age group to use Chiang's method. Wald's approximation for the variance of Mx needs at least 5 deaths for each age group (see Brown, Cai & DasGupta 2001).")
    }
    ax      <- rep(0.5, length(Mx))
    ax[1]   <- 0.1
    qx      <- Mx/(1 + (1 - ax) * Mx)
    qx[n]   <- 1
    qx[qx > 1] <- 1
    px      <- 1 - qx
    lx      <- c(1,cumprod(px)[-n])
    dx      <- -diff(lx)
    dx[n]   <- lx[n]
    
    Lx      <- c(lx[2:n] + ax[1:(n -  1)] * dx[1:(n - 1)], Mx[n] / ax[n])
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] <- 0
    Tx      <- rev(cumsum(rev(Lx)))
    ex      <- Tx/lx
    varqx   <- (Mx * (1 - ax * Mx)) / (Nx * (1 + (1 - ax) * Mx) ^3) # Chiang 2
    yx      <- lx^2 * ((1 - ax) + c(shift(ex,1)[-n],0)) ^2 * varqx
    cumyx   <- rev(cumsum(rev(yx)))
    varex   <- cumyx / (lx^2)
    seex    <- sqrt(varex)
    lower   <- ex[1] + qnorm((1 - level) / 2) * seex[1]
    upper   <- ex[1] - qnorm((1 - level) / 2) * seex[1]
  }
  
  cat("Loss of life expectancy due to the hump (years):",loss)
  cat("\n")
  cat("Half-confidence interval of fitted life expectancy at birth:",round((e0-lower)/2,3))
  cat("\n")
  if(loss > (e0-lower)/2){cat("The hump is statistically significant")}else{
    cat("The hump is not statistically significant")
  }
  cat("\n")
  cat("Significance level:",level)
  
  return(list(e0 = e0, upper = upper, lower = lower))
  
}

predict.morthump <- function(object, x = NULL, which = "fitted", tol = NULL,...){
  
  fit <- object
  
  model <- attr(fit, "model")
  
  if(is.null(x)){x <- fit$data$x}
  
  if(model == "hp"){
    
    if(which == "fitted"){y <- hp(x = x, par = coef(fit))}
    if(which == "hump" ){
      par <- coef(fit)
      par[-(4:6)] <- 0
      y <- hp(x = x, par = par)}
    
  }
  
  if(model == "hps"){
    
    if(which == "fitted"){y <- hps(x = x, par = coef(fit))}
    if(which == "hump" ){y <- rep(0, length(x))}
    
  }
  
  if(model == "hpk"){
    
    if(which == "fitted"){y <- hpk(x = x, par = coef(fit))}
    if(which == "hump" ){
      par <- coef(fit)
      par[-c(4:6,9)] <- 0
      y <- hpk(x = x, par = par)}
    
  }
  
  if(model == "sse"){
    
    decim <- Vectorize(function(x){min(which(x*10^(0:20) == floor(x*10^(0:20)))) - 1})
    
    if(all(x == round(x))){k <- 0}else{k <- max(decim(x))} # number of decimals required for the given input ages
    
    if(is.null(tol) == FALSE){k <- min(c(k,tol))} # limited if required
    
    n <- length(x)
    
    if(n == 1){x <- x + c(- 10^(-k), 0 , 10^(-k))} # for some reasons, the XX %*% coef step won't take a single value
    
    outer1 <- min(fit$data$x) # limits need to be identical to the original model to preserve the base
    outer2 <- max(fit$data$x)
    inner1 <- min(x)
    inner2 <- max(x)
    
    x1 <- unique(c(outer1, seq(inner1, inner2, 10^-k), outer2))
    
    xl <- min(x1)
    xr <- max(x1)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    
    deg2 <- fit$deg2
    deg3 <- fit$deg3
    
    B1 <- cbind(1,1/x1)
    B2 <- MortSmooth_bbase(x1, xmin, xmax, deg2, 3)
    B3 <- MortSmooth_bbase(x1, xmin, xmax, deg3, 3)
    XX <- list(X1=B1, X2=B2, X3=B3)
    
    XX <- lapply(XX, function(mat){mat[which(round(x1,k) %in% round(x,k)),]}) # only keep the ages we are interested in
    
    coef <- fit$coef
    nx <- length(XX)
    nc <- unlist(lapply(XX, ncol))
    ind2 <- cumsum(nc)
    ind1 <- c(1, ind2[1:(nx-1)]+1)
    ind <- NULL
    for(i in 1:nx){
      ind[[i]] <- ind1[i]:ind2[i]
    }
    
    etas <- NULL
    for(i in 1:nx){
      etas[[i]] <- XX[[i]] %*% coef[ind[[i]]]
    }
    eta1.hat <- etas[[1]]
    eta2.hat <- etas[[2]]
    eta3.hat <- etas[[3]]
    
    mhat <- lapply(etas, exp)
    
    if(n == 1){mhat <- lapply(mhat, function(mat){mat[2,]})}
    
    if(which == "fitted"){
      y <- rowSums(do.call(cbind,mhat))}
    if(which == "hump"){
      y <- mhat[[3]]
    }
    
  }
  
  return(as.vector(y))
  
}

LT <-
  function(Nx=NULL, Dx=NULL, Mx = Dx/Nx, ages = 0:(length(Mx)-1), axmethod = "midpoint", sex = "female",
           mxsmooth = TRUE, axsmooth = TRUE, radix = 1, verbose = TRUE){
    # the verbose function:
    Verb <- function(v, x){
      if (v) {
        cat(paste0(x, "\n"))
      }
    }
    # first a series of checks, messages and imputations to make sure the given Dx&Nx *or* Mx values are viable
    if (length(Mx) == 0){
      # two checks that will stop the function
      # 1) in absence of Mx, need both Dx and Nx
      if (is.null(Nx) | is.null(Dx)){
        Verb(verbose,"you're missing the key arguments")
        stop("either specify Mx directly, or else specify both Nx and Dx.")
      }
      # 2) both Nx and Dx must be of equal length
      stopifnot(length(Nx) != length(Dx))
      
      # safety to avoid zeros in the denominator. will make Mx of 1 in those cases
      if (any(Nx == 0)) {
        Verb(verbose, "there was at least one 0 in your Nx vector, imputed with corresponding Dx.")
        Nx[Nx == 0] <- Dx[Nx == 0]
      }
      Mx <- Dx / Nx
    }
    
    # by this point we should have an Mx vector with no NAs: no more need for Nx,Dx
    # we want to be able to accept 0s...
    
    
    # N is just used for counting, to save space
    N                   <- length(Mx) #
    Widths              <- diff(ages)
    Widths              <- c(Widths, Widths[N - 1])
    # define character for Age in formatted lifetable
    if (all(Widths == 1)){
      Age <- as.character(ages)
      Age[N] <- paste0(Age[N],"+")
    } else {
      Age <- c(0,paste(ages[2:(N - 1)], ages[2:(N - 1)] + Widths - 1,sep="-"),paste0(ages[N],"+"))
    }
    
    ages.mids.pre   	<- ages + Widths / 2
    ages.mids.pre[1] 	<- .1
    
    if(!is.null(Nx) & !is.null(Dx) & mxsmooth){
      # Giancarlo's package. I recommend either supplying an already-smoothed Mx vector (for complete control)
      # or else supplying Dx and Nx and setting mxsmooth to TRUE.
      # define some reasonable weights to block out certain cells if necessary
      w <- ifelse(Dx/Nx == 0 | is.null(Dx/Nx) | is.na(Dx/Nx) | log(Nx) < 0, 0, 1)
      fitBIC 			<- Mort1Dsmooth(x = ages.mids.pre, y = Dx, offset = log(Nx), w = w)
      Mx[2:N] 		<- (fitted(fitBIC) / Nx)[2:N]
    }
    
    if(is.null(Nx) & is.null(Dx) & mxsmooth){
      Verb(verbose,"mxsmooth was specified as TRUE, but since Mx was supplied directly, \nthere are no implicit weights (Nx). Function used a loess smoother \nto smooth out the Mx, but please be wary.")
      span 			<- ifelse(N > 30, .15, .4)
      logMx 			<- log(Mx)
      logMx[is.infinite(logMx)] <- NA
      logMx[2:N] 		<- predict(loess(logMx ~ ages.mids.pre, span = span, control = loess.control(surface = "interpolate")), newdata = ages[2:N])
      Mx 				<- exp(logMx)
    }
    
    # Assume that the lifetable mx is equal to the central death rate.
    mx                  <- Mx
    # these later 2 imputations should not be needed, but are there just in case.
    mx[is.na(mx)] 		<- 0
    mx[is.infinite(mx)] <- 0
    
    # we don't want 0s in m(x) for calculating a(x), because a(x) (except midpoint) is iterative
    # and we'd erroneously bring down neighboring age's ax values with zeros.
    # for later calulations, the zeros are taken 'as-is'
    Ind0 <- NULL
    if (length(axmethod) == 1 & min(mx) == 0){
      Ind0        <- mx == 0
      
      Verb(verbose, paste("\n\n*there were some ages (", ages[Ind0],
                          ") with no mortality.\nValues are imputed for calculating a(x), but the zeros are kept for the rest of the lifetable.\n"))
      span        <- ifelse(N > 30,.15,.4)
      logMx       <- log(mx)
      logMx[is.infinite(logMx)] <- NA
      Imp         <- exp(predict(loess(logMx ~ ages.mids.pre, span = span, control = loess.control(surface = "interpolate")), newdata = ages))[Ind0]
      if (any(is.na(Imp))){
        Imp     <- exp(spline(ages.mids.pre, logMx, xout=ages.mids.pre)$y[Ind0])
      }
      
      mx[Ind0]    <- Imp
    }
    ax <- NULL
    if (length(axmethod) == 1){
      if (axmethod %in% c("keyfitz", "schoen", "midpoint", "preston")){
        # here, the ax iteration is externalized to axEstimate()
        if (mxsmooth){
          Verb(verbose,"\nif mxsmooth = TRUE, then we don't smooth ax\n")
          axsmooth <- FALSE
        }
        ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = axmethod, sex = sex, verbose = verbose)
      }
      if (axmethod == "keyfitz" & !all(Widths == 1)){
        Verb(verbose, "\nIt appears you have an abridged lifetable, but have specified the keyfitz method of ax estimation.\nBe aware that this method presumes equal age intervals, as Preston et. al. (2001)\n warn on page 45. Consider using a different method or else specifying your own ax vector.\n Function continued nonetheless.")
      }
    }
    
    if (is.numeric(axmethod) & length(axmethod) == N) {
      ax <- axmethod
    }
    # last default # sex = "male"
    if (is.null(ax)){
      ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = "midpoint", sex = sex, verbose = verbose)
      Verb(verbose, "axmethod must be specified either as 'schoen','keyfitz','midpoint'\nor as a numeric vector the same length as Nx.\nThis wasn't the case, so the function defaulted to the midpoint method.")
    }
    
    # if zeros were imputed for ax estimation, then we put them back for the rest of the calculations
    if (!is.null(Ind0)){
      mx[Ind0]        <- 0
    }
    
    # now start the lifetable columns:
    qx 			        <- (Widths * mx) / (1 + (Widths - ax) * mx)
    qx[N] 		        <- 1
    
    # can't have qx > 1, so we impute 1s where necessary: hopefully only at the penultimate, as the case may be
    qx[qx > 1] 	        <- 1
    px 			        <- 1 - qx
    
    lx                  <- c(radix, radix * cumprod(px))[1:N]
    dx 			        <- c(-diff(lx),lx[N])
    
    Lx 	                <- c(Widths[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)], lx[N] * ax[N])
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] 	    <- 0
    
    Tx 			        <- rev(cumsum(rev(Lx)))
    ex 			        <- Tx / lx
    ex[N]               <- ifelse(mx[N] == 0, ax[N], {1 / mx[N]})
    
    # Sx is the pertinent output for projection matrices
    Sx   	            <- c((Lx[2:N] / Widths[2:N]) / (Lx[1:(N-1)] / Widths[1:(N - 1)]), Tx[N] / Tx[(N - 1)])
    
    # two not-very satisfying, and possibly redundant, substitutions:
    Sx[Lx == 0]   	    <- 0
    Sx[is.na(Sx)] 	    <- 0
    
    LT <- data.frame(cbind("Age" = Age, "ages" = ages, "mx" = round(mx, 5), "ax" = round(ax, digits = 2),
                           "qx" = round(qx, 5), "px" = round(px, 5), "lx" = round(lx, 5),
                           "dx" = round(dx, 5), "Lx" = round(Lx, 5), "Tx" = round(Tx, 5), "ex" = round(ex, 3)))
    # both LT as well as the individual pieces (not rounded) can be called
    out <- list(LT = LT,
                Age = Age,
                ages = ages,
                mx = mx,
                ax = ax,
                qx = qx,
                lx = lx,
                dx = dx,
                Lx = Lx,
                Tx = Tx,
                ex = ex,
                Sx = Sx,
                Widths = Widths)
    invisible(out)
  }

plot.codhump <- function(x, which, horiz = FALSE, ...){
  
  fit <- x
  
  if(which %in% 1:5 == FALSE){warning("Choose one of the five plots available.")}
  
  if(which == 1){
    x <- fit$mxc$x
    mxcd <- fit$mxc[,-c(1:2,ncol(fit$mxc))]
    mxwo <- fit$mxc$mxwo
    mx <- fit$mxc$mx
    col <- fit$par$col
    
    matplot(x, log(mxcd), type = "l", lty = 1, ylim = c(min(log(mxwo)),max(log(mx))), xlab = "age", ylab = "mx (log scale)", las = 1, col = col, axes = FALSE,...)
    lines(x, log(mxwo), col = "grey", lwd = 2)
    lines(x, log(mx), lwd = 2)
    legend(title = "Cause-deleted", x = "topleft", legend = names(fit$par$typ), col = fit$par$col, lty = 1, cex = 0.7)
    legend(title = "Remaining", x = "bottomright", legend = c("All-cause","Other causes"), lwd = 2, col = c(1,8), cex = 0.7)
    axis(1)
    axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
    box()
    title("Age-specific death rates")
  }
  
  
  if(which == 2){
    
    typ <- fit$par$typ
    fits <- fit$start
    x <- fit$mxc$x
    mxc <- fit$mxc$mx - fit$mxc[,-c(1:2,ncol(fit$mxc))]
    ncod <- length(typ)
    
    par(mfrow = disp(ncod, horiz = horiz))
    
    for(i in 1:length(fits)){
      
      plot(x, log(fit$mxc[,i+2]), col = 8, pch = 16, xlab = "age", ylab = "mx (log scale)", las = 1, axes = FALSE,...)
      lines(x, fits[[i]]$eta, col=2, lwd=2)
      lines(x, fits[[i]]$eta1, col=3)
      lines(x, fits[[i]]$eta2, col=4)
      #legend("top", legend = c("obs","total","hump","sen."), col=c(8,2:4),lty=c(-1,1,1,1), pch=c(16,-1,-1,-1),ncol=2)
      title(names(typ)[i])
      axis(1)
      axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
      box()
      
    }
    par(mfrow = c(1,1))
  }
  
  if(which == 3){
    
    dif <- fit$dif
    neg <- fit$neg
    
    par(mfrow = c(1,2))
    # speed of convergence
    plot(1:sum(dif!=0), dif[1:sum(dif!=0)], ylim = c(0,0.1), las = 1, xlab = "iteration", ylab = "max relative change of eta",...)
    abline(h = 1e-3, col = 2)
    title("Speed of convergence")
    
    
    # elimination of negative contributions
    neg <- neg[neg != -Inf & neg != Inf]
    plot(1:sum(neg!=0), - neg[1:sum(neg!=0)] * 100, ylim = c(0,-min(neg*100)), las = 1, xlab = "iteration", ylab = "% of negative values in the contributions",...)
    abline(h = 0, col = 2)
    title("Elimination of negative values")
    
    par(mfrow = c(1,1))
  }
  
  if(which == 4){
    
    # comparison of all-cause and cause-deleted fits
    typ <- fit$par$typ 
    x <- fit$mxc$x
    col <- fit$par$col
    ncod <- length(typ)
    fits <- fit$fits
    mxcd <- fit$mxc[,3:ncol(fit$mxc)]
    
    par(mfrow = disp(ncod, horiz = horiz))
    for(i in 1:ncod){
      fit.i <- fits[[i]]
      
      plot(x, log(fit$mxc[,2]), pch = 16, col = 1, xlab = "age", ylab = "mx (log scale)", ylim = c(-10,-1), axes = FALSE,...)
      points(x, log(mxcd[,i]), pch = 16, col = col[i])
      lines(x, fit$all$eta)
      lines(x, fit$all$eta1, lty = 2)
      lines(x ,fit$all$eta2, lty = 3)
      lines(x, fit.i$eta, col = col[i])
      lines(x, fit.i$eta1, col = col[i], lty = 2)
      lines(x, fit.i$eta2, col = col[i], lty = 3)
      polygon(x = c(x,rev(x)), y = c(fit$all$eta1,rev(fit.i$eta1)), col = paste(col[i],60,sep=""), border = NA)
      axis(1)
      axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
      box()
      title(names(typ)[i])
    }
    
    par(mfrow=c(1,1))
  }
  
  if(which == 5){
    
    typ <- fit$par$typ
    x <- fit$mxc$x
    col <- fit$par$col
    decomp <- fit$decomp
    decomp[decomp < 0] <- 0
    
    # cause- and age-specific contributions to the hump
    plot(x,fit$all$gamma1 / sum(fit$all$gamma1), xlim = c(10,60), type = "n", xlab = "age", ylab = "pdf", las = 1)
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(decomp[,1])) / sum(fit$all$gamma1),col = col[1], border = NA)
    if(ncol(decomp) > 1){
      for(i in 2:length(typ)){
        polygon(c(x,rev(x)),c(rowSums(as.matrix(decomp[,1:(i-1)])),rev(rowSums(as.matrix(decomp[,1:i])))) / sum(fit$all$gamma1),col = col[i], border = NA)
      }}
    lines(x,fit$all$gamma1 / sum(fit$all$gamma1),lwd = 2)
    
    #     axis(1,at = modeT, label = expression(M^T), lwd = 2, font = 2, line = 0.7)
    #     axis(1,at = modeA, label = expression(M^A), lwd = 3, col = col[1], font = 2, line = 0.7)
    #     axis(1,at=modeC,label=expression(M^C),lwd=3,col=col[3],font=2,line = 0.7)
    #     segments(x0=modeT,y0=0,x1=modeT,y1=max(fit$all$gamma1 / sum(fit$all$gamma1)),lty=2,lwd=2)
    #     segments(x0=modeA,y0=0,x1=modeA,y1=(fit$all$gamma1 / sum(fit$all$gamma1))[modeA-9],lty=2,col=col[1],lwd=2)
    #     segments(x0=modeC,y0=0,x1=modeC,y1=max(fit$decomp$C / sum(fit$all$gamma1)),lty=2,col=col[3],lwd=2)
    #     legend("topright", box.lty = 0, x.intersp = 0.3,
    #           legend = c("","area",expression(Delta^{kappa}),expression(bar(x)^kappa),expression(M^kappa),expression(sigma^kappa),
    #                      "A",paste(round(SA*100),"%"),paste(round(e10HA,2),"yr"),round(meanA,1),modeA,round(sdA,1),"C",paste(round(SC*100),"%"),paste(round(e10HC,2),"yr"),round(meanC,1),modeC,round(sdC,1)))
    
    legend("topright", fill = col, legend = names(typ))
    title("Contributions to the hump")
    
    #   negs <- pos <- decomp
    #   negs[negs > 0] <- 0
    #   pos[pos < 0] <- 0
    #   bp <- barplot(t(pos),space=0,col=col,border="grey",las=1,ylim=c(min(rowSums(negs)),max(rowSums(pos))))
    #   par(new = T)
    #   bp <- barplot(t(negs),space=0,col=col,border="grey",las=1,ylim=c(min(rowSums(negs)),max(rowSums(pos))),axes=F)
    #   #bp <- barplot(t(matrix(rep(mu1T,ncod),nrow=n) - do.call(cbind,muhat1)),space=0,col=col,border="grey",las=1)
    #   lines(bp-0.5, mu1T, lwd=2, type="s")
    #   axis(1,labels = x, at = x - min(x) + 0.5)
    #   title("Contributions to the hump: constrained")
    #   legend("topright",fill=col,legend=lab,border="grey")
    
  }
  
}

axEstimate <-
  function(Mx, n, axsmooth = TRUE, method = "keyfitz", sex, verbose){
    if (! method %in% c("keyfitz", "schoen", "midpoint", "preston")){
      stop("ax method specified is not valid, must be 'keyfitz','schoen','midpoint' or 'preston'")
    }
    if (method == "keyfitz"){
      ax <- axKeyfitz(Mx, n, axsmooth)
    }
    if (method == "schoen"){
      ax <- axSchoen(Mx, n, axsmooth)
    }
    if (method == "midpoint"){
      ax <- axMidpoint(Mx, n)
    }
    if (method == "preston"){
      if (missing(sex)) {
        cat(verbose,"\nWarning: You didn't specify 'sex', assumed 'female'!\n")
        sex <- "female"
      }
      ax <- axPreston(Mx, n, axsmooth, sex)
    }
    # use Andreev-Kingkade a0 formula.
    ax[1]  <- AKm02a0(Mx[1],sex)
    return(ax)
  }

axSchoen <-
  function(Mx, n, axsmooth = TRUE){
    N       <- length(Mx)
    if (axsmooth){
      ages    <- cumsum(n) - n
      span    <- ifelse(N > 30, .15, .4)
      Mx      <- log(Mx)
      Mx[2:N] <- predict(loess(Mx ~ ages,
                               span = span,
                               control = loess.control(surface = "interpolate")
      ),
      newdata = ages[2:N]
      )
      Mx      <- exp(Mx)
    }
    ax      <- ux   <- wx   <- lx   <- vector(length = N)
    lx[1]   <- 1
    for (i in 2:N){
      lx[i] <- lx[i - 1] * exp(-n[i - 1] * Mx[i - 1])
    }
    dx      <- -diff(lx)
    dx      <- c(dx, 1 - sum(dx))
    for (i in 2:(N - 1)){
      ux[i] <- ((n[i] ^ 2) / 240) * (Mx[i + 1] + 38 * Mx[i] + Mx[i - 1])
      wx[i] <- ((n[i] ^ 2) / 240) * (14 * Mx[i + 1] + 72 * Mx[i] - 6 * Mx[i - 1])
      ax[i] <- (ux[i] * lx[i] + wx[i] * lx[i + 1]) / dx[i]
    }
    ax[1]   <- .07 + 1.7 * Mx[1]
    ########### linear assumption for last ax ####################
    ########### using 3 info points           ####################
    coefs   <- lm(ax[(N - 3):(N - 1)] ~ c((N - 3):(N - 1)))$coef
    ax[N]   <- coefs[1] + N * coefs[2]
    return(ax)
  }

axPreston <-
  function(Mx, n, axsmooth = TRUE, sex = "female"){
    N <- length(Mx)
    if (axsmooth){
      ages    <- cumsum(n) - n
      span    <- ifelse(N > 30, .15, .4)
      Mx      <- log(Mx)
      Mx[2:N] <- predict(loess(Mx ~ ages,
                               span = span,
                               control = loess.control(surface = "interpolate")
      ),
      newdata = ages[2:N]
      )
      Mx      <- exp(Mx)
    }
    ax      <- n + (1 / Mx) - n / (1 - exp(-n * Mx))
    
    # Coale-Demeny rules of thumb:
    ax[1] <- ifelse(Mx[1] >= .107,
                    ifelse(sex == "female", .35, .33),
                    ifelse(sex == "female", {.053 + 2.8 * Mx[1]}, {.045 + 2.684 * Mx[1]})
    )
    if (n[2] == 4){
      ax[2] <- ifelse(Mx[1] >= .107,
                      ifelse(sex == "female", 1.352, 1.361),
                      ifelse(sex == "female", {1.522 - 1.581 * Mx[1]}, {1.651 - 2.816 * Mx[1]})
      )
    }
    # my own adaptation for a1-a8 in the case of single years: otherwise they're *very* close to .5
    if (length(unique(n[1:10])) == 1){
      if (unique(n[1:10]) == 1){
        int         <- ax[9] - ax[1]
        jumps       <- c(100, 30, 10, 8, 6, 4, 2)
        jumps       <- jumps / sum(jumps)
        for (i in 1:7){
          ax[(i + 1)]     <- ax[i] + jumps[i] * int
        }
      }
    }
    return(ax)
  }

disp <- function(n, horiz = FALSE){
  
  r <- c <- floor(sqrt(n))
  
  if(horiz == FALSE){
    while(r*c < n){r <- r + 1}}else{
      while(r*c < n){c <- c + 1}}
  
  return(c(r,c))
}

axMidpoint <-
  function(Mx, n){
    ax      <- .5 * n
    ax[1]   <- .07 + 1.7 * Mx[1] # I guess this is Keyfitz + Flieger?
    return(ax)
  }

axKeyfitz <-
  function(Mx, n, axsmooth = TRUE){
    # iterative ax-dx process decribed on page 44-45 of
    # Preston et al, Demography: Measuring and Modelling Population Processes. Blackwell Publishing, 2001
    N <- length(Mx)
    if (axsmooth){
      ages        <- cumsum(n) - n
      span        <- ifelse(N > 30, .15, .4)
      Mx          <- log(Mx)
      Mx[2:N]     <- predict(loess(Mx ~ ages,
                                   span = span,
                                   control = loess.control(surface = "interpolate")
      ),
      newdata = ages[2:N]
      )
      Mx          <- exp(Mx)
    }
    axit        <- .5 * n
    axit[1] <- .07 + 1.7 * Mx[1]
    for (i in 1:7){
      qx              <- (n * Mx) / (1 + (n - axit) * Mx)
      qx[length(Mx)]  <- 1
      px              <- 1 - qx
      lx              <- 1 # typically radix would go here, but it makes no difference since values don't pass on.
      for (i in 2:length(Mx))	{
        lx[i]   <- lx[i - 1] * px[i - 1]
      }
      dx          <- -diff(lx)
      for (i in 2:(length(Mx) - 1)){
        axit[i] <- (-(n[i - 1] / 24) * dx[i - 1] + (n[i] / 2) * dx[i] + (n[i + 1] / 24) * dx[i + 1]) / dx[i]
      }
      
      # this is just my own way of finishing off the ax's, not sooo creative,
      # but it doesn't usually make a difference
      axit[N - 1] <- axit[N - 2] - (axit[N - 3] - axit[N - 2]) * 1.5
      axit[N]     <- axit[N - 1] - (axit[N - 2] - axit[N - 1]) * 1.5
      # it assumes continued senescence at the final ages:
      axit[N - 1] <- axit[N - 2] - (axit[N - 3] - axit[N - 2]) * 1.5
      axit[N]     <- axit[N - 1] - (axit[N - 2] - axit[N - 1]) * 1.5
    }
    axit[1]     <- .07 + 1.7 * Mx[1]
    return(axit)
  }

codgroup <- function(data, x.range = 10:35, k = "ASW"){
  
  x <- data$x
  mxc <- as.data.frame(data$mxc)
  lab <- as.data.frame(data$lab)
  names(mxc) <- data$lab$short
  names(mxc)[is.na(names(mxc))] <- LETTERS[1:sum(is.na(names(mxc)))]
  
  rxc <- as.data.frame(apply(mxc,2,diff))
  names(rxc) <- names(mxc)
  d <- dist(x = t(rxc[which(x > min(x.range) & x < max(x.range)),]),method = "euclidian")
  cluster <- hclust(d)
  
  opt <- as.clustrange(object = cluster, diss = d, ncluster = ncol(rxc))
  if(!is.numeric(k)){
    k <- summary(opt)[which(rownames(summary(opt)) == k),1]
  }
  
  pca <- prcomp(t(rxc[which(x > min(x.range) & x < max(x.range)),]))
  
  gr <- cutree(cluster, k = k)
  
  typ <- list()
  for(i in 1:k){typ[[i]] <- which(gr == i)}
  names(typ) <- LETTERS[1:k]
  typ <- typ[-which.max(lapply(typ,length))]
  
  out <- list(cluster = cluster, groups = gr, k = k, typ = typ, data = data, x.range = x.range)
  
  class(out) <- "codgroup"
  
  return(out)
  
}

summary.codhump <- function(object, ...){
  
  fit <- object
  
  # loss of life expectancy : all-cause
  mhat <- fit$all$gamma1 + fit$all$gamma2
  mhyp <- mhat - fit$all$gamma1
  widths <- diff(fit$mxc$x) ; widths <- c(widths[1],widths)
  if(max(widths) == 1){ages <- fit$mxc$x}else{ages <- floor(fit$mxc$x/widths)*widths}
  e0hat <- LT(Mx = mhat, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hyp <- LT(Mx = mhyp, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hump  <- round(e0hyp - e0hat,2)
  
  # loss of life expectancy : cause-specific
  losses <- rep(NA, length(fit$par$typ))
  for(i in 1:length(fit$par$typ)){
    mhump <- fit$all$gamma1 - fit$fits[[i]]$gamma1
    mhyp <- mhat - mhump
    e0hyp <- LT(Mx = mhyp, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
    losses[i] <- round(e0hyp - e0hat, 2)
  }
  losses <- c(e0hump,losses)
  names(losses) <- c("all",names(fit$par$typ))
  
  # descriptive statistics
  x <- fit$mxc$x
  mh <- as.numeric(fit$all$gamma1)
  mode <- x[which.max(mh)]
  mean <- mh %*% x / sum(mh)
  sd <- sqrt(mh %*% (x - c(mean))^2 / sum(mh))
  modes <- means <- sds <- rep(NA, length(fit$par$typ))
  for(i in 1:length(modes)){
    modes[i] <- x[which.max(fit$decomp[,i])]
    means[i] <- fit$decomp[,i] %*% x / sum(fit$decomp[,i])
    suppressWarnings(sds[i]   <- sqrt(fit$decomp[,i] %*% (x - c(means[i]))^2 / sum(fit$decomp[,i])))
  }
  modes <- c(mode,modes)
  means <- c(mean,means)
  sds <- c(sd, sds)
  modes <- round(modes, 2)
  means <- round(means, 2)
  sds <- round(sds, 2)
  names(sds) <- names(means) <- names(modes) <- c("all",names(fit$par$typ))
  
  out <- list(loss = losses, mode = modes, mean = means, sd = sds)
  
  class(out) <- "summary.codhump"
  attr(out, "extra") <- fit
  
  invisible(out)
  
}

#' @rdname summary.codhump
#' @method print summary.codhump
#' @export

print.summary.codhump <- function(x, ...){
  
  fit <- attr(x, "extra")
  data <- fit$par$data
  share.abs <- function(c){round(sum(data$dxc[data$x >= 10 & data$x < 35,c]) / sum(data$dxc[data$x >= 10 & data$x < 35,]) * 100, 1)}
  
  cat("[[>>]] Information about the model used to fit the cause- and age-specific death rates")
  cat("\nNumber of causes contributing to the hump:",length(fit$par$typ))
  cat("\nList of causes contributing to the hump:",names(fit$par$typ))
  cat("\nShare of observed deaths between 10 and 34 (%):")
  cat("\n",paste0(names(fit$par$typ)," = ",unlist(lapply(fit$par$typ, share.abs)),collapse = ", "))
  cat("\nAge range used for the analysis:",min(fit$par$x.range),"-",max(fit$par$x.range))
  cat("\nNumber of iterations:",length(fit$dif),"( max =",fit$par$maxit,")")
  
  if(length(fit$dif) == fit$par$maxit){lab <- "Reached maximum number of iterations"}
  if(diff(rev(fit$dif))[1] < 1e-3){lab <- "Reached a stable solution"}
  if(fit$neg[length(fit$neg)] < fit$neg[length(fit$neg)-1]){lab <- "Could not approach the constraints any closer"}
  cat("\nReason for stopping the optimization:",lab)
  if(any(x$loss[-1] / x$loss[1] < 5e-2) | max(fit$neg) < -0.01){
    cat("\nConsider changing parameters to improve fit (see documentation from codhump)")
    cat(paste("\n>> Small contribution of :",names(fit$par$typ)[x$loss[-1] / x$loss[1] < 5e-2]))
    cat(paste("\n>> Negative contributions amount to",-round(max(fit$neg*100),1),"% of total hump"))
  }
  
  
  cat("\n\n")
  
  cat("[[>>]] Descriptive statistics about the hump")
  
  cat("\nIntensity")
  cat("\nLoss of life expectancy due to the hump: ", x$loss[1], "year")
  cat("\nLoss of life expectancy by cause (years):",paste0(paste0(names(fit$par$typ)," = "),x$loss[-1], collapse = ", "))
  cat("\nLoss of life expectancy by cause (%):",paste0(paste0(names(fit$par$typ)," = "), round(x$loss[-1] / x$loss[1] * 100, 1), collapse = ", "))
  
  cat("\nInaccuracy compared to overall hump:",round((sum(x$loss[-1]) - x$loss[1]) / x$loss[1] * 100, 1),"%")
  cat("\n")
  
  cat("\nCentrality")
  cat("\nMode of overall hump:", round(x$mode[1],2))
  cat("\nMode by cause:",paste0(paste0(names(fit$par$typ)," = "),x$mode[-1], collapse = ", "))
  cat("\nMean of overall hump:", round(x$mean[1],2))
  cat("\nMean by cause:",paste0(paste0(names(fit$par$typ)," = "),round(x$mean[-1],2), collapse = ", "))
  cat("\n")
  
  cat("\nDispersion")
  cat("\nStandard deviation of overall hump:", round(x$sd[1],2))
  cat("\nStandard deviation by cause:",paste0(paste0(names(fit$par$typ)," = "),round(x$sd[-1],2), collapse = ", "))
  
  
  
  invisible(x)
  
}

unabridge <- function(dxc, nx, inter, plot = FALSE){
  
  fun <- function(v, .inter){
    
    x <- c(.inter[1,1],.inter[,2])
    y <- c(0,cumsum(v))
    
    if(v[1] > 0.5 * sum(v)){
      
      xx <- min(x):(max(x)-1)
      yy <- rep(v / apply(.inter,1,diff), times = apply(.inter,1,diff))
      ycum <- c(0,cumsum(yy))
      
    }else{
      
      spl <- smooth.spline(x = x, y = y, all.knots = TRUE, spar = -1)
      ycum <- predict(spl, x = min(x):max(x))$y
      
      yy <- diff(ycum)
      xx <- min(x):(max(x)-1)
      
    }
    
    if(plot == TRUE){
      par(mfrow = c(1,2))
      plot(x, y, type = "s", las = 1, xlab = "age", ylab = "cdf")
      lines(min(x):max(x), ycum, lty = 3, col = 2, type = "s")
      
      plot(.inter[,1], diff(y), type = "s", xlab = "age", ylab = "pdf", las = 1)
      lines(xx, yy * rep(apply(.inter,1,diff), times = apply(.inter,1,diff)),
            type = "s", col = 2, lty = 3)
      par(mfrow = c(1,1))
    }
    
    return(yy)
    
  }
  
  dxc.unabr <- apply(dxc[-1,],2,fun,inter[-1,])
  
  dxc.unabr[dxc.unabr < 0] <- 0
  
  dxc.unabr <- rbind(dxc[1,],dxc.unabr)
  
  dxc.unabr <- as.data.frame(dxc.unabr)
  
  names(dxc.unabr) <- names(dxc)
  
  nx.unabr <- fun(nx, inter)
  nx.unabr[nx.unabr < 0] <- 0
  
  x <- inter[1,1]:(max(inter[,2])-1)
  
  return(list(dxc = dxc.unabr, x = x, nx = nx.unabr))
  
}

hps <- function(x, par){
  
  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  G <- par[5]
  H <- par[6]
  
  y <- A ^((x + B) ^ C) + D  + (G * (H ^ x))/(1 + G * (H ^ x))
  
  return(y)
}

xtp <- function(path, max, country, year, sex, cmx){
  
  if(sex %in% c("males","females") == F){warning("Extrapolation of non-extinct cohorts is currently only available for males and females.")}
  
  per <- read.demogdata(file = paste(path,country,"STATS","Mx_1x1.txt",sep="/"), popfile = paste(path,country,"STATS","Exposures_1x1.txt",sep="/"), type = "mortality", label = country)
  
  per$age <- 0:max
  per$rate$male <- per$rate$male[1:(max+1),]
  per$rate$female <- per$rate$female[1:(max+1),]
  per$pop$male <- per$pop$male[1:(max+1),]
  per$pop$female <- per$pop$female[1:(max+1),]
  
  if(any(colSums(is.na(per$rate$male)) > 10)){
    
    problem <- which(colSums(is.na(per$rate$male)) > 10)
    per$rate$male[,problem] <- per$rate$male[,min(problem)-1]
    per$rate$female[,problem] <- per$rate$female[,min(problem)-1]
    per$pop$male[,problem] <- per$pop$male[,min(problem)-1]
    per$pop$female[,problem] <- per$pop$female[,min(problem)-1]
    
  }
  
  end <- rev(per$year)[1]
  nonext <- (end - max):(end - 30)
  
  fun <- function(x){sum(is.na(x)) != 0 & is.na(x[1]) == F}
  av <- by(data = cmx$Female, INDICES = cmx$Year, FUN = fun)
  gens <- unique(cmx$Year)[av == T]
  
  cMx <- as.data.frame(matrix(nrow = max+1, ncol = length(gens)), row.names = 0:max)
  names(cMx) <- gens
  
  # projection
  s <- ifelse(sex == "males", "male", "female")
  fdm <- fdm(per, series = s)
  fc <- demography::forecast.fdm(fdm, h = max-30+1)
  
  if(sex == "males"){cmx.fc <- diag(fc$rate$male[(end-year+1):(max+1),])}
  if(sex == "females"){cmx.fc <- diag(fc$rate$female[(end-year+1):(max+1),])}
  
  cmx <- cmx[cmx$Year == year,]
  cmx[(end-year+1):(max+1),ifelse(sex == "males",4,3)] <- cmx.fc
  
  return(cmx)
  
}

AKq02a0 <- function(q0, sex = "m"){
  sex <- rep(sex, length(q0))
  ifelse(sex == "m",
         ifelse(q0 < .0226, {0.1493 - 2.0367 * q0},
                ifelse(q0 < 0.0785, {0.0244 + 3.4994 * q0},.2991)),
         ifelse(q0 < 0.0170, {0.1490 - 2.0867 * q0},
                ifelse(q0 < 0.0658, {0.0438 + 4.1075 * q0}, 0.3141))
  )
}

AKm02q0 <- function(m0, constant, slope){
  -1 / slope / m0 * (-m0 +  (m0 * constant) - 0.1e1 + sqrt(((m0 ^ 2) - 2 * constant * (m0 ^ 2) + 2 * m0 + (constant ^ 2) * (m0 ^ 2) - 2 *  (m0 * constant) + 1 - 4 * slope * (m0 ^ 2)))) / 2
}

AKm02a0 <- function(m0,sex="male"){
  sex <- rep(sex,length(m0))
  ifelse(sex == "male",
         ifelse(m0 < 0.02306737, 0.1493 - 2.0367 * AKm02q0(m0, 0.1493, -2.0367),
                ifelse(m0 < 0.0830706, 0.0244 + 3.4994 * AKm02q0(m0, 0.0244, 3.4994), .2991)),
         ifelse(m0 < 0.01725977, 0.1490 - 2.0867 * AKm02q0(m0, 0.1490, -2.0867),
                ifelse(m0 < 0.06919348, 0.0438 + 4.1075 * AKm02q0(m0, 0.0438, 4.1075), 0.3141))
  )
}

cod.data <- function(k, Mxc, exp){
  
  if(!(all(k < 0) | all(k > 0))){break ; warning("Cod must be all included or all excluded")}
  
  if(length(k) == 1){if(k > 0){m <- Mxc[,k]}else{m <- rowSums(Mxc[,k])}}else{
    m <- rowSums(Mxc[,k])}
  
  if(!is.null(dim(m))){warning("Error")}
  
  d <- m * exp
  
  return(data.frame(m,d,exp))
}

hps.fit <- function(data, method = "port", w = 1/data$m, start = NULL){
  
  # according to Heligman and Pollard (1980), weights should be 1/q
  # according to Brillinger (1986), weights should be 1/d
  # according to Heligman and Pollard (1980), the response variable should be q or q/(1-q), but it is defined here as m
  
  
  if(data$x[1] == 0){data$x[1] <- 1e-5}
  
  form <- as.formula(m ~ A ^ (x + B) ^ C + D + G * H ^ x / (1 + G * (H ^ x)))
  
  if(is.null(start)){
    
    start <- list(A = 0.001, B = 0.5, C = 0.11, D = 0.0015, G = 3e-5, H = 1.105)
    
    lower <- c(0.0001, 0.000001, 0.0001, 0, 0.0000001, 0.5)
    
    upper <- c(0.1, 0.5, 1, 0.01, 0.01, 1.5)
    
  }else{
    lower <- start$lower
    upper <- start$upper
    start <- start$start
  }
  
  if(method %in% c("port","lm","gnm","bayes") == FALSE){warning("The method must be one of 'port, 'lm' or 'gnm'.")}
  
  if(method == "port"){
    
    fit <- nls(formula = form, data = data, start = start, lower = lower, upper = upper,
               algorithm = "port", weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "lm"){
    
    fit <- nlsLM(formula = form, data = data, start = start, lower = lower, upper = upper,
                 weights = w, control = list(maxiter=1000))
    fit$coef <- coef(fit)
  }
  
  if(method == "gnm"){
    
    # see Currie 2014 and Debon et al. 2005
    
    warning("Sorry, the estimation of the Siler-Heligman-Pollard model as a generalized linear model is not implemented yet.")
    
  }
  
  if(method == "bayes"){
    
    # fun <- function(m){rnorm(n = 8e3, mean = m, sd = m/10)}
    
    # prior <- do.call(cbind,lapply(start,fun))
    
    # fit <- hp.bm.imis(prior = prior, nrisk = data$exp, ndeath = data$d, K = 10)
    
    warning("Sorry, the estimation of the Siler-Heligman-Pollard model using Bayesian statistics is not implemented yet.")
    
  }
  
  if(data$x[1] < 1){data$x[1] <- 0}
  
  fit$data <- data
  fit$method <- method
  fit$w <- w
  fit$type <- "parametric"
  
  
  return(fit)
  
}