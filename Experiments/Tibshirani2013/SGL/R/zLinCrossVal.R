linCrossVal<-
function(data, index, nfold = 10, nlam = 20, min.frac = 0.05, alpha = 0.95, thresh = 0.0001, maxit = 10000, gamma = 0.8, verbose = TRUE, step = 1, reset = 10){

  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)

   ## Setting up group lasso stuff ##
  
  ord <- order(index)
  index <- index[ord]
  X <- X[,ord]
  unOrd <- match(1:length(ord),ord)

  ## Coming up with other C++ info ##

  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0,(num.groups+1))
  for(i in 1:num.groups){
    range.group.ind[i] <- min(which(index == groups[i])) - 1
  }
  range.group.ind[num.groups+1] <- ncol(X)

  group.length <- diff(range.group.ind)
  beta.naught <- rep(0,ncol(X))
  beta <- beta.naught

  ## Done with group stuff ##

  ## finding the path

MainSol <- oneDim(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, gamma = gamma, step = step, reset = reset, alpha = alpha)

  lambdas <- MainSol$lambdas

  lldiff <- rep(0, nlam)
  lldiffFold <- matrix(0, nrow = nlam, ncol = nfold)
  
  size <- ceiling(nrow(data$x)/nfold)
  ind <- sample(1:nrow(data$x), rep = FALSE)
  for(i in 1:nfold){
    if(i < nfold){
      ind.out <- ind[((i-1)*size+1):(i*size)]
      ind.in <- ind[-(((i-1)*size+1):(i*size))]
  }
    if(i == nfold){
      ind.out <- ind[((i-1)*size+1):n]
      ind.in <- ind[-(((i-1)*size+1):n)]
    }
  
    new.data <- list(x = data$x[ind.in,], y = data$y[ind.in])

    new.sol <- oneDim(new.data, index, thresh = thresh, inner.iter = maxit, lambdas = lambdas, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, gamma = gamma, step = step, reset = reset, alpha = alpha)

    for(k in 1:nlam){
    
	lldiffFold[k,i] <- sum((y[ind.out] - X[ind.out,] %*% new.sol$beta[ ,k])^2) / 2

        lldiff[k] <- lldiff[k] + sum((y[ind.out] - X[ind.out,] %*% new.sol$beta[ ,k])^2) / 2
      }
    if(verbose == TRUE){
  write(paste("*** NFOLD ", i, "***"),"")
  }
  }
  lldiffSD <- apply(lldiffFold,1,sd)*sqrt(nfold)
  obj <- list(lambdas = lambdas, lldiff = lldiff,llSD = lldiffSD, fit = MainSol)
  class(obj)="cv.SGL"
  return(obj)
}
