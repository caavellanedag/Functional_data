#
# Autor(es): Camilo Avellaneda
# Mantenimiento: Camilo Avellaneda
# Fecha creación: 22/04/2021
#==============================================

if(!("input" %in% list.files())){dir.create("input")}
if(!("output" %in% list.files())){dir.create("output")}
if(!("src" %in% list.files())){dir.create("src")}
if(!("docs" %in% list.files())){dir.create("docs")}

require(pacman)
p_load(tidyverse,openxlsx,here,haven,data.table,fda,fda.usc,roahd,fdaoutlier,caret)
p_load(maptools,sgeostat,scatterplot3d,car,fields,gstat,geoR)

options(scipen=999)
#lambda_10<-10^10
#y<-df[,.y]
#x<-range
#y = (x-20)^2+200*rnorm(1000)0

# Tabla de datos a trabajar -----------------------------------------------
data(tecator)

absorp_2 <- absorp %>%
  as.data.table() %>% t() %>% 
  data.table::melt(value.name="Absorbencia") %>% as.data.table()
absorp_2 <- absorp_2[,Medicion := as.numeric(str_replace_all(Var1,"V",""))]

## Se pegan los grupos de acuerdo a las categor?as de grasa
colnames(endpoints)<-paste0("V",1:3)
endpoints<-endpoints %>% as.data.table()
endpoints[,c("Grupo","Llave"):= list(ifelse(V2>20,"Alto (Grasa>20)","Bajo (Grasa<=20)"),1:nrow(endpoints))]
absorp_3 <- merge(absorp_2,endpoints[,c("Llave","V2","Grupo")],by.x="Var2",by.y="Llave",all.x=TRUE)


absorp_3 <- absorp_3[,Curva:=Var2]

#Gr?fico descriptivo de la realizaci?n
descriptive_1 <- absorp_2 %>% ggplot()+
  geom_line(aes(Medicion,Absorbencia,group=Var2),color="midnightblue")+
  guides(color=F)+theme_bw()+xlab("Medición")+ylab("Absorbancia")

#Gr?fico descriptivo de la realizaci?n por grupo
descriptive_2 <- absorp_3 %>% 
  ggplot()+
  geom_line(aes(Medicion,Absorbencia,group=Var2),color="midnightblue")+
  guides(color=F)+theme_bw()+xlab("Medición")+ylab("Absorbancia")+
  facet_grid(.~Grupo)



SelectionNumBasis <- function(x, Y_mat, Num_basis){
spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis = Num_basis)
object_fd <- smooth.basis(y=Y_mat, fdParobj=spline_basis)
predict_fd <- eval.fd(x, fdobj=object_fd$fd)
SSE <- sum((predict_fd-Y_mat)^2)
CME <- SSE/(length(x)-Num_basis)
X <- eval.basis(x, spline_basis)
Var_estimator <- sum(diag(X%*% solve(t(X) %*% X) %*% t(X)))*CME
return(data.frame(Num_basis=Num_basis,Bias_squared=SSE,Var_estimator=Var_estimator,MSE=Var_estimator+SSE))
}

Select_lambda <- function(x,Y_mat,lambda){
  lambda_10 <- 10^lambda
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=lambda_10)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  #mean((eval.fd(1:100, fdobj=absorb_fd$fd) - t(absorp)[,.y])^2)
  data.frame(lambda = lambda,gcv = sum(object_fd$gcv))
}


df_num_basis <- map_dfr(4:50,~SelectionNumBasis(x = 1:100, Y_mat = t(absorp),Num_basis = .x))
df_num_basis %>% melt(id.vars="Num_basis",variable.name="Variable",value.name="Value") %>% 
  ggplot()+geom_line(aes(x=Num_basis,Value,linetype=Variable,color=Variable),size=1.1)+
  scale_color_brewer(palette="Set1")+theme_bw()+ylab("Cuadrado medio del error total")+
  xlab("Número de funciones base")+ylim(0,5)

df_gcv <- map_dfr(seq(-5,5,by=0.01),~Select_lambda(x= x , Y_mat= Y_mat, lambda=.x)) 
lambda_optimo <- df_gcv[which.min(df_gcv$gcv),"lambda"]

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


Matrix_bajos <- t(absorp[as.matrix(endpoints[str_detect(Grupo,"Bajo"),"Llave"]),])
Matrix_altos <- t(absorp[as.matrix(endpoints[str_detect(Grupo,"Bajo"),"Llave"]),])
Bajos_fitted <- smooth_fda(x = 1:100, Y_mat = Matrix_bajos, lambda = lambda_optimo)
Altos_fitted <- smooth_fda(x = 1:100, Y_mat = Matrix_altos, lambda = lambda_optimo)

fD_bajo <- fData( 1:100, t(Bajos_fitted))
boxplot_bajo <- fbplot(fD)

fD_alto <- fData( 1:100, t(Altos_fitted))
boxplot_alto <- fbplot(fD)

Matrix_alto_without <- Matrix_altos[,- boxplot_alto$ID_outliers]
Matrix_bajo_without <- Matrix_bajos[,- boxplot_bajo$ID_outliers]


# Calcula la estadística de prueba
Test_change_point <- function(x, Y_mat, k, n_pca, lambda,...){
  N <- ncol(Y_mat)
  Y_mat_1 <- Y_mat[,1:k]
  Y_mat_2 <- Y_mat[,(k+1):ncol(Y_mat)]
  fitted_1 <- fitted_fda( x = x, Y_mat = Y_mat_1, lambda = lambda)
  fitted_2 <- fitted_fda( x = , Y_mat = Y_mat_2, lambda = lambda)
  fda_all <- fitted_fda( x = x, Y_mat = Y_mat, lambda = lambda)
  pca_overall <- pca.fd(fda_all$fd, nharm = n_pca)
  lambda <- pca_overall$values
  T_test <- sum((rowSums(inprod(pca_overall$harmonics, fitted_1$fd))-
  k/N*rowSums(inprod(pca_overall$harmonics, fitted_2$fd)))^2/lambda[1:n_pca])
  return(T_test)
}

# Calcula n_times veces la estadística de prueba de cambio de punto
p_value_boots_change_k <- function(x, Y_mat, n_pca, lambda, n_times, k = k){
test_stat_or <- map_dbl(1,~Test_change_point(x, Y_mat = Y_mat,
                                       lambda =  lambda, n_pca = n_pca, k = k))
test_stat_boots <- map_dbl(1:n_times,~Test_change_point(x, 
                        Y_mat = Y_mat[,sample(1:ncol(Y_mat), ncol(Y_mat) ,replace=T)],
                        lambda =  lambda, n_pca = n_pca, k = k))
mean(test_stat_or < test_stat_boots)
}

# Compila los resultados para diferentes valores en el número de funciones en cada grupo 
p_value_boots_change_all <- function(x, Y_mat, n_pca, lambda, n_times){
  map_dbl(seq(10,ncol(Y_mat)-10,by=5),
      ~p_value_boots_change_k(x, Y_mat, n_pca, lambda, n_times, k = .x))
}

# Se revisa si hay algún punto donde los p valores de p_value_boots_change_all donde se rechace la prueba  
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

Check_differences(x = 1:100, Y_mat = Matrix_bajo_without, n_pca =1, lambda = lambda_optimo, n_times=20)



Independence_Q <- function(x, Y_mat, lambda, n_pca){
  pca <- pca.fd(fitted_fda(x = x, Y_mat = Y_mat, lambda = lambda)$fd, nharm = n_pca)
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
  map_dbl(1:(ncol(Y_mat)-1),~Find_Q_h(.x)) %>% sum()
}

p_value_boots <- function(x, Y_mat, lambda, n_pca, n_times, FUN){
test_stat_or <- Independence_Q(x = x, Y_mat = Y_mat, lambda = lambda, n_pca = n_pca)
test_stat_boots <- map_dbl(1:n_times,~FUN(x,
                                          Y_mat = Y_mat[,sample(1:ncol(Y_mat), ncol(Y_mat) ,replace=T)],
                                          lambda =  lambda, n_pca = n_pca))
mean(test_stat_or < test_stat_boots)
}

p_value_boots(x = 1:100, Y_mat = Matrix_bajo_without, lambda = lambda_optimo, n_pca = 2, n_times = 20, FUN = Independence_Q)



# Primer prueba de Hipótesis - Prueba adecuada ----------------------------
# Se va a probar si la función promedio poblacional es constante e igual al promedio
# general del grupo con alto contenido de grasa



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



Test_one_mean_fda(x = 1:100, Y_mat = Matrix_bajo_without, lambda = lambda_optimo,
                  n_pca = 1, mean = mean(Matrix_bajo_without), n_times = 60)



## Prueba con respecto a dos poblaciones

Y_mat_1 =  Matrix_bajo_without
Y_mat_2 =  Matrix_alto_without
n_pca = 1
lambda = lambda_optimo
n_times <- 10
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

P_value_test_two_mean(x = 1:100, Y_mat_1 = Matrix_alto_without,
                      Y_mat_2 = Matrix_bajo_without, lambda = lambda_optimo,
                      n_pca = 1, n_times = 50)




Y_mat_1 =  Matrix_bajo_without
Y_mat_2 =  Matrix_alto_without
n_pca = 1
lambda = lambda_optimo
n_times <- 10


Y_mat <- cbind(Y_mat_1,Y_mat_2)
pca_overall <- pca.fd(fitted_fda(x = 1:100, Y_mat = Y_mat, lambda = lambda_optimo)$fd,
                      nharm = n_pca)

Cov_operators_test_stat <- function(x , Y_mat_1, Y_mat_2, pca_overall, lambda){
N <- ncol(Y_mat_1)
M <- ncol(Y_mat_2)
theta <- N/(N+M)
fitted_1 <- eval.fd(evalarg= x,
                    fitted_fda(x = x, Y_mat = Y_mat_1, lambda = lambda)$fd)
fitted_2 <- eval.fd(evalarg= x,
                    fitted_fda(x = x, Y_mat = Y_mat_2, lambda = lambda)$fd)
epsilon_x <- inprod(fitted_fda(x = x, Y_mat = Y_mat_1, lambda = lambda)$fd,
pca_overall$harmonics)

epsilon_y <- inprod(fitted_fda(x = x, Y_mat = Y_mat_2, lambda = lambda)$fd,
       pca_overall$harmonics)
lambda_x <- 1/N * map_dbl(1:1,~sum(epsilon_x[,.x]^2))
lambda_y <- 1/M * map_dbl(1:1,~sum(epsilon_y[,.x]^2))


C_x <- map_dfc(1:nrow(epsilon_x),function(.x){
epsilon_x[.x,1] * fitted_1[,.x]
}) %>% as.matrix() %>% rowMeans()

C_y <- map_dfc(1:nrow(epsilon_y),function(.x){
  epsilon_y[.x,1] * fitted_2[,.x]
}) %>% as.matrix() %>% rowMeans()


Numerator <- inprod(fitted_fda(x = x, Y_mat = C_x - C_y, lambda = lambda)$fd , pca_overall$harmonics)^2
T <- (N+M)/2 * theta * (1-theta) * (Numerator)/((theta * lambda_y + (1-theta)*lambda_x)*(theta * lambda_x + (1-theta)*lambda_y))
return(T)
}





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


Results_boots <- map_dfr(1:n_times,function(n){
  Simulated_matrix <- map_dfc(1:(N+M),function(m){
    if(runif(1) > N/(N+M)){
      epsilon_star <- Mean_2 + Res_2[,sample(1:M,size=1)]
    }else{
      epsilon_star <- Mean_1 + Res_1[,sample(1:N,size=1)]
    }
    return(epsilon_star)
  }) 
  
  
  Mean_2_sim <- eval.fd(evalarg = x, 
                        mean.fd(fitted_fda(x = x, Y_mat = as.matrix(Simulated_matrix[,(N+1):(N+M)]),
                                           lambda = lambda)$fd))
  Mean_1_sim <- eval.fd(evalarg = x, 
                        mean.fd(fitted_fda(x = x, Y_mat = as.matrix(Simulated_matrix[,1:N]),                                     
                                           lambda = lambda)$fd))
  
  Test_stat_boots <- Cov_operators_test_stat(x = 1:100 , Y_mat_1 = as.matrix(Simulated_matrix[,1:N]),
                          Y_mat_2 = as.matrix(Simulated_matrix[,(N+1):(N+M)]), 
                          pca_overall = pca_overall, 
                          lambda = lambda_optimo)
  
  
  return(data.frame(Test_stat_boots = Test_stat_boots))
})


Test_original <- Cov_operators_test_stat(x = 1:100 , Y_mat_1 = Matrix_bajo_without,
                                         Y_mat_2 = Matrix_alto_without, 
                                         pca_overall = pca_overall, 
                                         lambda = lambda_optimo)

p_value <- pchisq(Test_original, n_pca*(n_pca+1)/2, lower.tail = FALSE)
p_value

mean(Results_boots$Test_stat_boots > as.double(Test_original))
