#
# Autor(es): Camilo Avellaneda
# Mantenimiento: Camilo Avellaneda
# Fecha creación: 04/05/2021
#==============================================
getwd()
if(!("input" %in% list.files())){dir.create("input")}
if(!("output" %in% list.files())){dir.create("output")}
if(!("src" %in% list.files())){dir.create("src")}
if(!("docs" %in% list.files())){dir.create("docs")}
if(!("reports" %in% list.files())){dir.create("reports")}


source("Z:/Universidad Nacional/PhD/Functional_data/Taller 4/src/Funciones.R")

require(pacman)
p_load(tidyverse,openxlsx,here,haven,data.table,fda,fda.usc,roahd,fdaoutlier,caret)
p_load(maptools,sgeostat,scatterplot3d,car,fields,gstat,geoR)

options(scipen=999)

data(tecator)
endpoints<-endpoints %>% as.data.table()
endpoints[,c("Grupo","Llave"):= list(ifelse(V2>20,"Alto (Grasa>20)","Bajo (Grasa<=20)"),1:nrow(endpoints))]

k_optimo <- map_dbl(1:nrow(absorp),function(y){
  x <- 1:100
  k<-k.check(gam(absorp[y,]~s(x)))
  return(k[1,1])
}) %>% mean() %>% round()

Y <- endpoints$V2
fit.fpcr = pfr(Y ~ fpc(absorp,k=50))
# without a penalty, use k= 3 for the first beta, and 15 for the second
# one can verify the efficacy of these choices by
# looking at the aic
fit.lin = pfr(Y ~ lf(absorp, bs = "ps", k = k_optimo, fx = TRUE))
# "ps" stands for "penalized splines", fx= TRUE means no penalty is used
fit.pfr = pfr(Y ~ lf(absorp, bs = "ps", k = k_optimo))
# if sp is not specified, data driven smoothing is used
grid <- 1:100

coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr)$value,
                   Basis = coef(fit.lin)$value,
                   Penalized = coef(fit.pfr)$value)
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")


plot_1 <- ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                    = Method),width=12,height=6) + geom_path(size=1) + theme_bw()
plot_2 <- ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                    = Method),width=12,height=6) + geom_path(size=1) + theme_bw()+
                    ylim(-2*10^3,2*10^3)


x_axis <- 1:215
FPCR_predict <- predict(fit.fpcr)
Basis_predict <- predict(fit.lin)
Penalized_predict <- predict(fit.pfr)
True_values <- Y

df_to_plot <- data.frame(x_axis, FPCR_predict, Basis_predict, Penalized_predict, True_values)

plot_3 <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  ggplot()+geom_point(aes(x_axis,True_values),color="red")+
  geom_point(aes(x_axis,Estimado),color="blue")+
  facet_wrap(~Tipo)


Table_MSE <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  mutate(Squared_dif = (True_values-Estimado)^2) %>% group_by(Tipo) %>%
  summarize(MSE = mean(Squared_dif)) %>% ungroup()
  

Table_R_squared <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  mutate(Squared_dif = (True_values-Estimado)^2) %>% group_by(Tipo) %>%
  summarize(SSE = sum(Squared_dif)) %>% ungroup() %>% 
  mutate(SST = sum((Y-mean(Y))^2)) %>% mutate(R2 = 1-SSE/SST)





df_gcv <- map_dfr(seq(-5,5,by=0.01),~Select_lambda(x= 1:100 , Y_mat= base::t(absorp), lambda=.x)) 
lambda_optimo <- df_gcv[which.min(df_gcv$gcv),"lambda"]

absorp_fda <- fitted_fda(x = 1:100, Y_mat = t(absorp), lambda = lambda_optimo )  

deriv_fda <- deriv.fd(absorp_fda$fd)
deriv_matrix <- t(eval.fd(evalarg = 1:100, deriv_fda))



fit.fpcr_2 = pfr(Y ~ fpc(absorp,k=50) + fpc(deriv_matrix,k=50))
# without a penalty, use k= 3 for the first beta, and 15 for the second
# one can verify the efficacy of these choices by
# looking at the aic
fit.lin_2 = pfr(Y ~ lf(absorp, bs = "bs", k = k_optimo, fx = TRUE) +
                lf(deriv_matrix, bs = "bs", k = k_optimo, fx = TRUE))
# "ps" stands for "penalized splines", fx= TRUE means no penalty is used
fit.pfr_2 = pfr(Y ~ lf(absorp, bs = "bs", k = k_optimo) + 
                lf(deriv_matrix, bs = "bs", k = k_optimo))
# if sp is not specified, data driven smoothing is used



coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr_2)$value,
                   Basis = coef(fit.lin_2)$value,
                   Penalized = coef(fit.pfr_2)$value)
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")


plot_4 <- ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                              = Method),width=12,height=6) + geom_path(size=1) + theme_bw()
plot_5 <- ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                              = Method),width=12,height=6) + geom_path(size=1) + theme_bw()+
  ylim(-2*10^3,2*10^3)


x_axis <- 1:215
FPCR_predict <- predict(fit.fpcr_2)
Basis_predict <- predict(fit.lin_2)
Penalized_predict <- predict(fit.pfr_2)


df_to_plot <- data.frame(x_axis, FPCR_predict, Basis_predict, Penalized_predict, True_values)

plot_6 <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  ggplot()+geom_point(aes(x_axis,True_values),color="red")+
  geom_point(aes(x_axis,Estimado),color="blue")+
  facet_wrap(~Tipo)


Table_MSE <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  mutate(Squared_dif = (True_values-Estimado)^2) %>% group_by(Tipo) %>%
  summarize(MSE = mean(Squared_dif)) %>% ungroup()


Table_R_squared <- df_to_plot %>%
  data.table::melt(variable.name="Tipo",value.name="Estimado",id.vars=c("x_axis","True_values")) %>% 
  mutate(Squared_dif = (True_values-Estimado)^2) %>% group_by(Tipo) %>%
  summarize(SSE = sum(Squared_dif)) %>% ungroup() %>% 
  mutate(SST = sum((Y-mean(Y))^2)) %>% mutate(R2 = 1-SSE/SST)


