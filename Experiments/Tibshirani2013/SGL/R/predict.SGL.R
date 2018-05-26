predict.SGL = function(x,newX,lam){
  cvobj = x

  X <- newX

  if(!is.null(x$X.transform)){
    X <- t(t(newX) - x$X.transform$X.means)
    X <- t(t(X) / x$X.transform$X.scale)
  }

  intercept <- 0

  if(!is.null(x$intercept)){
    intercept <- x$intercept[lam]
  }

  eta <- X %*% x$beta[,lam] + intercept


  if(x$type == "linear"){
    y.pred <- eta
  }

  if(x$type == "logit"){
    y.pred = exp(eta)/(1+exp(eta))
  }

 if(x$type == "cox"){
    y.pred = exp(eta)
  }

  return(y.pred)
}
