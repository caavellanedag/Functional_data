Check_differences <- function(x, Y_mat, n_pca, lambda, n_times){
  p_values <- p_value_boots_change_all(x = x,
                                       Y_mat = Y_mat, n_pca = n_pca,
                                       lambda = lambda, n_times = n_times)
  K <- seq(10,ncol(Y_mat)-10,by=5)
  i=1
  x <- p_values[i]
  while(x > 0.05 & i != length(K)){
    i=i+1
    x <- p_values[i]
  }
  
  if(i != length(p_values)){
    Diferencia <- K[i]
  }else{
    Diferencia <- ncol(Y_mat)
  }
  return(list(K = K, p_values = p_values, Diferencia = Diferencia))
}

Find_C_h <- function(h){
  C_h <- matrix(0, nrow = n_pca, ncol = n_pca)
  for(k in 1:n_pca){
    for(m in 1:n_pca){
      C_h[k,m] = pca$scores[1:(nrow(pca$scores)-h),k] %*% pca$scores[(h+1):(nrow(pca$scores)),m]
    }
  }
  return(C_h)
}


Find_Q_h <- function(h){
  C_h <- matrix(0, nrow = n_pca, ncol = n_pca)
  for(k in 1:n_pca){
    for(m in 1:n_pca){
      C_h[k,m] = pca$scores[1:(nrow(pca$scores)-h),k] %*% pca$scores[(h+1):(nrow(pca$scores)),m]
    }
  }
  Q_h <- sum(diag(base::t(C_h) %*% solve(diag(pca$values[1:n_pca])) %*% C_h %*% solve(diag(pca$values[1:n_pca]))))
  return(Q_h)
}

fitted_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  return(object_fd)
}

Independence_Q <- function(x, Y_mat, lambda, n_pca){
  pca <- pca.fd(fitted_fda(x = x, Y_mat = Y_mat, lambda = lambda)$fd, nharm = n_pca)
  Find_Q_h <- function(h){
    C_h <- matrix(0, nrow = n_pca, ncol = n_pca)
    for(k in 1:n_pca){
      for(m in 1:n_pca){
        C_h[k,m] = pca$scores[1:(nrow(pca$scores)-h),k] %*% pca$scores[(h+1):(nrow(pca$scores)),m]
      }
    }
    if(n_pca == 1){
      Q_h <- sum(diag(base::t(C_h) %*% 1/pca$values[n_pca] %*% C_h %*% 1/pca$values[1:n_pca]))    
    }else{
      Q_h <- sum(diag(base::t(C_h) %*% solve(diag(pca$values[1:n_pca])) %*% C_h %*% solve(diag(pca$values[1:n_pca]))))  
    }
    
    return(Q_h)
  }
  map_dbl(1:(ncol(Y_mat)-1),~Find_Q_h(.x)) %>% sum()
}

p_value_boots <- function(x, Y_mat, lambda, n_pca, n_times, FUN){
  test_stat_or <- Independence_Q(x = x, Y_mat = Y_mat, lambda = lambda, n_pca = n_pca)
  test_stat_boots <- map_dbl(1:n_times,~FUN(x,
                                            Y_mat = Y_mat[,sample(1:ncol(Y_mat), ncol(Y_mat) ,replace=T)],
                                            lambda =  lambda, n_pca = n_pca))
  mean(test_stat_or < test_stat_boots)
}

p_value_boots_change <- function(x, Y_mat, n_pca, lambda, n_times, k = k){
  test_stat_or <- map_dbl(1,~Test_change_point(x, Y_mat = Y_mat,
                                               lambda =  lambda, n_pca = n_pca, k = k))
  test_stat_boots <- map_dbl(1:n_times,~Test_change_point(x, 
                                                          Y_mat = Y_mat[,sample(1:ncol(Y_mat), ncol(Y_mat) ,replace=T)],
                                                          lambda =  lambda, n_pca = n_pca, k = k))
  mean(test_stat_or < test_stat_boots)
}



p_value_boots_change_all <- function(x, Y_mat, n_pca, lambda, n_times){
  map_dbl(seq(10,ncol(Y_mat)-10,by=5),
          ~p_value_boots_change_k(x = x, Y_mat = Y_mat, n_pca =n_pca,
                                  lambda = lambda, n_times = n_times, k =.x))
}

p_value_boots_change_k <- function(x, Y_mat, n_pca, lambda, n_times, k = k){
  test_stat_or <- map_dbl(1,~Test_change_point(x, Y_mat = Y_mat,
                                               lambda =  lambda, n_pca = n_pca, k = k))
  test_stat_boots <- map_dbl(1:n_times,
                             ~Test_change_point(x,
                                                Y_mat = Y_mat[,sample(1:ncol(Y_mat), ncol(Y_mat) ,replace=T)],
                                                lambda =  lambda, n_pca = n_pca, k = k))
  mean(test_stat_or < test_stat_boots)
}


p_value_boots_one_mean <- function(x, Y_mat, n_pca, lambda, n_times){
  test_stat_or <- Test_one_mean_fda(x = x, 
                                    Y_mat = Y_mat, 
                                    lambda = lambda, 
                                    n_pca = n_pca,  mean = mean(Y_mat))
  
  n <- ncol(Y_mat)
  Matrix_alto_test_1 <- fitted_fda(x = x, Y_mat = Y_mat, lambda = lambda_optimo)
  Boots_matrix <- eval.fd(evalarg = 1:100, fdobj = center.fd(Matrix_alto_test_1$fd))+mean(Y_mat)
  
  test_stat_boots <- map_dfr(1:n_times,~Test_one_mean_fda(x, 
                                                          Y_mat = Boots_matrix[,sample(1:ncol(Boots_matrix),ncol(Boots_matrix),replace=TRUE)],
                                                          lambda =  lambda, n_pca = n_pca, mean = mean(Y_mat)))
  #mean(test_stat_or > test_stat_boots)
  return(data.frame(p_value_norm = mean(test_stat_boots$T_norm > test_stat_or$T_norm),
                    p_value_fpca = mean(test_stat_boots$T_fpca > test_stat_or$T_fpca)))
}


Select_lambda <- function(x,Y_mat,lambda){
  lambda_10 <- 10^lambda
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=lambda_10)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  #mean((eval.fd(1:100, fdobj=absorb_fd$fd) - t(absorp)[,.y])^2)
  data.frame(lambda = lambda,gcv = sum(object_fd$gcv))
}


SelectionNumBasis <- function(x, Y_mat, Num_basis){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis = Num_basis)
  absorb_fd <- smooth.basis(y=Y_mat, fdParobj=spline_basis)
  predict_fd <- eval.fd(x, fdobj=absorb_fd$fd)
  SSE <- sum((predict_fd-Y_mat)^2)
  CME <- SSE/(length(x)-Num_basis)
  X <- eval.basis(x, spline_basis)
  Var_estimator <- sum(diag(X%*% solve(base::t(X) %*% X) %*% base::t(X)))*CME
  return(data.frame(Num_basis=Num_basis,Bias_squared=SSE,Var_estimator=Var_estimator,MSE=Var_estimator+SSE))
}

smooth_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  predict_fd <- eval.fd(x, fdobj=object_fd$fd)  
  return(predict_fd)
}

Stats_boots <- function(){
  fitted_fda_boots <- fitted_fda(x = x, Y_mat = Boots_matrix[,sample(1:ncol(Boots_matrix),ncol(Boots_matrix),replace=TRUE)] - mean(Y_mat) ,
                                 lambda = lambda_optimo)
  muhat_boots <- mean.fd(fitted_fda_boots$fd)
  T_norm_boots <- n*sum(inprod(v,muhat_boots)^2)
  T_fpca_boots <- n*sum(inprod(v,muhat_boots)^2/lambda)
  return(data.frame(T_norm = T_norm_boots, T_fpca = T_fpca_boots))
}

Test_change_point <- function(x, Y_mat, k, n_pca, lambda,...){
  N <- ncol(Y_mat)
  Y_mat_1 <- Y_mat[,1:k]
  Y_mat_2 <- Y_mat[,(k+1):ncol(Y_mat)]
  fitted_1 <- fitted_fda( x = x, Y_mat = Y_mat_1, lambda = lambda)
  fitted_2 <- fitted_fda( x = x, Y_mat = Y_mat_2, lambda = lambda)
  fitted <- fitted_fda( x = x, Y_mat = Y_mat, lambda = lambda)
  fda_all <- fitted_fda( x = x, Y_mat = Y_mat, lambda = lambda)
  pca_overall <- pca.fd(fda_all$fd, nharm = n_pca)
  lambda <- pca_overall$values
  T_test <- sum((rowSums(inprod(pca_overall$harmonics, fitted_1$fd))-
                   k/N*rowSums(inprod(pca_overall$harmonics, fitted$fd)))^2/lambda[1:n_pca])
  return(T_test)
}


Test_one_mean_fda <- function(x, Y_mat, lambda, n_pca, mean, n_times){
  n <- ncol(Y_mat)
  Matrix_alto_test_1 <- fitted_fda(x = x, Y_mat = Y_mat - mean, lambda = lambda)
  pca_test_1 <- pca.fd(Matrix_alto_test_1$fd, nharm = n_pca,centerfns = FALSE)
  
  v <- pca_test_1$harmonics
  lambda <- pca_test_1$values[1:n_pca]
  muhat <- mean.fd(Matrix_alto_test_1$fd)
  
  T_norm_or <- n*sum(inprod(v,muhat)^2)
  T_fpca_or <- n*sum(inprod(v,muhat)^2/lambda)
  
  Boots_matrix <- eval.fd(evalarg = x, 
                          fdobj = center.fd(fitted_fda(x = x, Y_mat = Y_mat, lambda = lambda_optimo)$fd))+
    mean
  
  Stats_boots <- function(x, Y_mat, Boots_matrix, v, lambda, mean){
    fitted_fda_boots <- fitted_fda(x = x, Y_mat = Boots_matrix[,sample(1:ncol(Boots_matrix),ncol(Boots_matrix),replace=TRUE)] - mean ,
                                   lambda = lambda)
    muhat_boots <- mean.fd(fitted_fda_boots$fd)
    T_norm_boots <- n*sum(inprod(v,muhat_boots)^2)
    T_fpca_boots <- n*sum(inprod(v,muhat_boots)^2/lambda)
    return(data.frame(T_norm = T_norm_boots, T_fpca = T_fpca_boots))
  }
  
  final_stats_boots <- map_dfr(1:n_times,~Stats_boots(x, Y_mat, Boots_matrix, v, lambda, mean))
  
  
  return(data.frame(p_value_norm = mean(final_stats_boots$T_norm>T_norm_or),
                    p_value_fpca = mean(final_stats_boots$T_fpca>T_fpca_or)))
}



P_value_test_two_mean <- function(x, Y_mat_1, Y_mat_2, lambda, n_pca, n_times){
  Y_mat = cbind(Y_mat_1, Y_mat_2) 
  N <- ncol(Y_mat_1)
  M <- ncol(Y_mat_2)
  
  Cov_overall <- N/(N+M)*cov(base::t(Y_mat_1))+
    M/(N+M)*cov(base::t(Y_mat_2))
  
  Mean_2 <- eval.fd(evalarg = x, 
                    mean.fd(fitted_fda(x = x, Y_mat = Y_mat_2, lambda = lambda)$fd))
  Mean_1 <- eval.fd(evalarg = x, 
                    mean.fd(fitted_fda(x = x, Y_mat = Y_mat_1, lambda = lambda)$fd))
  Res_1 <- eval.fd(evalarg = x, 
                   fitted_fda(x = x, Y_mat = Y_mat_1, lambda = lambda)$fd) - 
    (Mean_1 %*% matrix(1,nrow = 1, ncol =N))
  
  Res_2 <- eval.fd(evalarg = x, 
                   fitted_fda(x = x, Y_mat = Y_mat_2, lambda = lambda)$fd) - 
    (Mean_2 %*% matrix(1,nrow = 1, ncol =M))
  
  
  Mean <- eval.fd(evalarg = x, 
                  mean.fd(fitted_fda(x = x, Y_mat = Y_mat, lambda = lambda)$fd))
  
  lambda_2 <- eigen(Cov_overall)$values[1:n_pca]
  v <- eigen(Cov_overall)$vectors[,1:n_pca]
  
  epsilon_or <- inprod(fitted_fda(x = x, Y_mat = Mean_1 - Mean_2, lambda = lambda)$fd,
                       fitted_fda(x = x, v, lambda = lambda)$fd)
  T_norm_or <- N*M/(N+M)*epsilon_or^2
  T_fpca_or <- sum(epsilon_or^2/lambda_2)
  
  Results_boots <- map_dfr(1:n_times,function(n){
    Simulated_matrix <- map_dfc(1:(N+M),function(m){
      if(runif(1) > N/(N+M)){
        epsilon_star <- Mean + Res_2[,sample(1:M,size=1)]
      }else{
        epsilon_star <- Mean + Res_1[,sample(1:N,size=1)]
      }
      return(epsilon_star)
    })
    
    
    Mean_2_sim <- eval.fd(evalarg = x, 
                          mean.fd(fitted_fda(x = x, Y_mat = as.matrix(Simulated_matrix[,(N+1):(N+M)]),
                                             lambda = lambda)$fd))
    Mean_1_sim <- eval.fd(evalarg = x, 
                          mean.fd(fitted_fda(x = x, Y_mat = as.matrix(Simulated_matrix[,1:N]),                                     
                                             lambda = lambda)$fd))
    epsilon_sim <- inprod(fitted_fda(x = x, Y_mat = Mean_1_sim - Mean_2_sim,
                                     lambda = lambda)$fd,
                          fitted_fda(x = x, v, lambda = lambda)$fd)
    T_norm_boots <- sum(N*M/(N+M)*epsilon_sim^2)
    T_fpca_boots <- sum(epsilon_sim^2/lambda_2)
    return(data.frame(T_norm_boots = T_norm_boots, T_fpca_boots = T_fpca_boots))
  })
  
  return(data.frame(p_value_norm = mean(Results_boots$T_norm_boots > as.double(T_norm_or)),
                    p_value_fpca = mean(Results_boots$T_fpca_boots > as.double(T_fpca_or))))
}





