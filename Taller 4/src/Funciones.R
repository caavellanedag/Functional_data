
smooth_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  predict_fd <- eval.fd(x, fdobj=object_fd$fd)  
  return(predict_fd)
}


fitted_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  return(object_fd)
}


Select_lambda <- function(x,Y_mat,lambda){
  lambda_10 <- 10^lambda
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=lambda_10)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  #mean((eval.fd(1:100, fdobj=absorb_fd$fd) - t(absorp)[,.y])^2)
  data.frame(lambda = lambda,gcv = sum(object_fd$gcv))
}

