library(fda)
library(tidyverse)
library(data.table)
library(caret)
library(RColorBrewer)
options(scipen=999)
date<-str_replace(str_replace_all(Sys.Date(),"-",""),"20","")
setwd("E:/Documentos Camilo/PhD/Functional data analysis")

SelectionNumBasis<-function(.x,.y){
  #lambda_10<-10^10
  y<-t(absorp)[,.y]
  #y = (x-20)^2+200*rnorm(1000)0
  x = 1:100
  # put noise to the true curve and generate noisy data
  spline_basis=create.bspline.basis(rangeval=c(1,100),nbasis=.x)
 # fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^(-1.75))
  #fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=0)
  absorb_fd=smooth.basis(y=y, fdParobj=spline_basis)
  Bias_2<-t(eval.fd(x, fdobj=absorb_fd$fd) - y) %*% (eval.fd(x, fdobj=absorb_fd$fd) - y) 
  CME<-Bias_2/(100-.x)
  X<-eval.basis(1:100, spline_basis)
  Var_estimator<-sum(diag(X%*% solve(t(X) %*% X) %*% t(X)))*CME
  data.frame(Num_function=.y,Num_basis=.x,Bias_squared=Bias_2,Var_estimator=Var_estimator,MSE=Var_estimator+Bias_2)
}

data(tecator)




df<-map2_dfr(rep(4:35,each=nrow(absorp)),rep(1:nrow(absorp),32),~SelectionNumBasis(.x,.y))
df<-data.table(df)

df[,.(MSE=sum(MSE),Bias_squared=sum(Bias_squared),
      Var_estimator=sum(Var_estimator)),
   by=.(Num_basis)] %>% melt(id.vars="Num_basis",variable.name="Variable",value.name="Value") %>% 
  ggplot()+geom_line(aes(x=Num_basis,Value,linetype=Variable,color=Variable),size=1.1)+
  scale_color_brewer(palette="Set1")+theme_bw()+ylab("Cuadrado medio del error total")+
  xlab("Número de funciones base")+ylim(0,5)
Summary_num_basis<-df[,.(MSE=sum(MSE),Bias_squared=sum(Bias_squared),
      Var_estimator=sum(Var_estimator)),
   by=.(Num_basis)]

Summary_num_basis[order(MSE),]
OptimalNumBasis<-34



absorp_2<-absorp %>%
  as.data.table() %>% t() %>%
  data.table::melt(value.name="Absorbencia") %>% as.data.table()
absorp_2<-absorp_2[,Medicion := as.numeric(str_replace_all(Var1,"V",""))]


#scale_color_manual(values=colorRampPalette(brewer.pal(8, "Set1"))(215))

colnames(endpoints)<-paste0("V",1:3)
endpoints<-endpoints %>% as.data.table()
endpoints[,c("Grupo","Llave"):= list(ifelse(V2>20,"Alto (Grasa>20)","Bajo (Grasa<=20)"),1:nrow(endpoints))]
absorp_3<-merge(absorp_2,endpoints[,c("Llave","V2","Grupo")],by.x="Var2",by.y="Llave",all.x=TRUE)


absorp_3<-absorp_3[,Curva:=Var2]




spline_basis=create.bspline.basis(rangeval=c(1,100),nbasis=34)
fdParobj<- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^(-1.75))
ggplot_fd<-function(df,iter){
  df.fd=smooth.basis(y=absorp_3[Curva==iter,Absorbencia], fdParobj=fdParobj)
  predict(df.fd,newdata=1:100)
}

predict_fd<-map_dfc(1:length(unique(absorp_3$Curva)),
                    ~ggplot_fd(df=absorp_3,iter=.x)) %>% as.data.table()
colnames(predict_fd)<-paste0("X",1:length(unique(absorp_3$Curva)))


predict_fd<-predict_fd[,c("labels"):=list(1:100)] %>% 
  melt(variable.name="Variable",value.name="Fitted",id.vars="labels")  
predict_fd<-predict_fd[,c("Variable"):=list(as.numeric(str_replace_all(Variable,"X","")))]
absorp_4<-merge(absorp_3,predict_fd,by.x=c("Curva","Medicion"),by.y=c("Variable","labels"),all.x=TRUE)





absorp_4_altos<-absorp_4[Grupo=="Alto (Grasa>20)",c("Curva","Medicion","Absorbencia")]
absorp_4_bajos<-absorp_4[Grupo=="Bajo (Grasa<=20)",c("Curva","Medicion","Absorbencia")]

absorp_4_altos<-dcast(absorp_4_altos,Medicion~Curva,value.var="Absorbencia")
absorp_4_altos<-absorp_4_altos[,-"Medicion"]

absorp_4_bajos<-dcast(absorp_4_bajos,Medicion~Curva,value.var="Absorbencia")
absorp_4_bajos<-absorp_4_bajos[,-"Medicion"]

fit_fd_altos=smooth.basis(y=as.matrix(absorp_4_altos), fdParobj=fdParobj)
fit_fd_bajos=smooth.basis(y=as.matrix(absorp_4_bajos), fdParobj=fdParobj)

mean_fd_altos=mean.fd(fit_fd_altos$fd)
trimmed_mean_fd_altos=mean.fd(fit_fd_altos$fd,trim=0.1)
std_fd_altos=std.fd(fit_fd_altos$fd)

mean_fd_bajos=mean.fd(fit_fd_bajos$fd)
trimmed_mean_fd_bajos=mean.fd(fit_fd_bajos$fd,trim=0.1)
std_fd_bajos=std.fd(fit_fd_bajos$fd)

summary_stats<-data.frame(
  rbind(
data.frame(Medicion=1:100,Mean=predict(mean_fd_bajos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Bajo (Grasa<=20)"),
data.frame(Medicion=1:100,Mean=predict(mean_fd_altos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Alto (Grasa>20)")))

predict(std_fd_altos,newdata=1:100)

ggplot(absorp_4)+
  geom_line(data=absorp_4 ,aes(x=Medicion,y=Fitted,group=Curva),color="darkorange",size=1,alpha=0.14)+
  geom_line(data=summary_stats,aes(x=Medicion,y=mean),size=2,color="darkred",linetype="twodash")+
  geom_line(data=summary_stats,aes(x=Medicion,y=Sd),size=2,color="green",linetype="twodash")+
  guides(color=FALSE)+theme_bw()+xlab("Medición")+ylab("Absorbencia (Curvas ajustadas)")+
  facet_grid(.~Grupo)

# mean_fd=predict(mean_pinch,newdata=1:nrow(pinch))
# sd_fd=predict(sd_pinch,newdata=1:nrow(pinch))
# summary_stats=data.frame(labels=1:nrow(pinch),mean_fd,sd_fd) %>% as.data.table() %>% 
#   melt(id.vars="labels",value.name="Value",variable.name="Stat")


var_fd_bajos<-var.fd(fit_fd_bajos$fd)
var_fd_altos<-var.fd(fit_fd_altos$fd)

grid=(1:100)
W.cov.mat=eval.bifd(grid, grid, var_fd_bajos)
persp(grid, grid, W.cov.mat, xlab="s",
      ylab="t", zlab="c(s,t)")
contour(grid, grid, W.cov.mat, lwd=2)

W.cov.mat=eval.bifd(grid, grid, var_fd_altos)
persp(grid, grid, W.cov.mat, xlab="s",
      ylab="t", zlab="c(s,t)")
contour(grid, grid, W.cov.mat, lwd=2)

GCV<-map_dfc(1:ncol(t(absorp)),function(.y){
  map_dbl(seq(-5,5,by=0.05),function(.x){
    lambda_10<-10^.x
    spline_basis=create.bspline.basis(rangeval=c(1,100),nbasis=34)
    fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=lambda_10)
    absorb_fd=smooth.basis(y=t(absorp)[,.y], fdParobj=fdParobj)
    #mean((eval.fd(1:100, fdobj=absorb_fd$fd) - t(absorp)[,.y])^2)
    absorb_fd$gcv
  })
}) 


df_gcv<-data.frame(log10_lambda=seq(-5,5,by=0.05),
lambda=10^log10_lambda,
GCV=rowMeans(GCV))
df_gcv[which.min(df_gcv$GCV),]



.y<-75
.x<-34
  #lambda_10<-10^10
  y<-t(absorp)[,.y]
  #y = (x-20)^2+200*rnorm(1000)0
  x = 1:100
  # put noise to the true curve and generate noisy data
  spline_basis=create.bspline.basis(rangeval=c(1,100),nbasis=.x)
  fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^(-1.75))
  #fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=0)
  absorb_fd=smooth.basis(y=y, fdParobj=fdParobj)
  # Bias_2<-t(eval.fd(x, fdobj=absorb_fd$fd) - y) %*% (eval.fd(x, fdobj=absorb_fd$fd) - y) 
  CME<-Bias_2/(100-.x)
  X<-eval.basis(1:100, spline_basis)
  Var_estimator<-diag(X%*% solve(t(X) %*% X) %*% t(X))*CME
  
  return(data.frame(Prediction=mean.fd(absorb_fd$fd)
                    Error=sqrt(Var_estimator),
                    LI=mean.fd(absorb_fd$fd)-200*sqrt(Var_estimator),
                    LS=mean.fd(absorb_fd$fd)+200*sqrt(Var_estimator)))
  #data.frame(Num_function=.y,Num_basis=.x,Bias_squared=Bias_2,Var_estimator=Var_estimator,MSE=Var_estimator+Bias_2)
}
selection_curve_df<-as.data.table(Selection_curve_confidence_band(.x=34,.y=75))
selection_curve_df<-selection_curve_df[,labels:=1:100]
ggplot()+geom_line(data=selection_curve_df,aes(labels,Prediction))+
  geom_line(data=selection_curve_df,aes(labels,LI),linetype="twodash")+


descriptive_1<-absorp_2 %>% ggplot()+geom_line(aes(Medicion,Absorbencia,group=Var2),color="midnightblue")+
  guides(color=F)+theme_bw()+xlab("Medición")+ylab("Absorbancia")

descriptive_2<-absorp_3 %>% 
  ggplot()+
  geom_line(aes(Medicion,Absorbencia,group=Var2),color="midnightblue")+
  guides(color=F)+theme_bw()+xlab("Medición")+ylab("Absorbancia")+
  facet_grid(.~Grupo)



tiff(paste0(date,"_Descriptive_1_.tiff"), units="in", width=7, height=7, res=200)
print({
  descriptive_1
})
dev.off()


tiff(paste0(date,"_Descriptive_2_.tiff"), units="in", width=7, height=7, res=200)
print({
  descriptive_2
})
dev.off()


tiff(paste0(date,"_Selection_num_basis_.tiff"), units="in", width=7, height=7, res=200)
print({
  df[,.(MSE=sum(MSE),Bias_squared=sum(Bias_squared),
      Var_estimator=sum(Var_estimator)),
   by=.(Num_basis)] %>% melt(id.vars="Num_basis",variable.name="Variable",value.name="Value") %>% 
  ggplot()+geom_line(aes(x=Num_basis,Value,linetype=Variable,color=Variable),size=1.1)+
  scale_color_brewer(palette="Set1")+theme_bw()+ylab("Cuadrado medio del error total")+
  xlab("Número de funciones base")+ylim(0,5)
})
dev.off()


tiff(paste0(date,"_Selection_lambda_.tiff"), units="in", width=7, height=7, res=200)
print({
df_gcv %>% ggplot()+geom_point(aes(x=log10_lambda,GCV),size=1)+theme_bw()+
  ylab("Validación cruzada generalizada")+xlab("Log10(lambda)")
})
dev.off()




absorp_selected_curve<-absorp_4[,c("Curva","Medicion","Absorbencia")] %>% dcast(Medicion~Curva,value.var="Absorbencia")
absorp_selected_curve<-absorp_selected_curve[,-"Medicion"]
fit_selected_curve=smooth.basis(y=as.matrix(absorp_selected_curve), fdParobj=fdParobj)
std_fd=sd.fd(fit_selected_curve$fd)

fdParobj = fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^(-1.75))
selected_curve_fd<-smooth.basis(y=as.matrix(absorp_3[Curva==75,"Absorbencia"]), fdParobj=fdParobj)
df_selected<-data.frame(fitted=predict(selected_curve_fd),se=predict(std_fd,newdata=1:100)) %>% as.data.table()
var(as.matrix(absorp_3[Curva==75,"Absorbencia"]))

tiff(paste0(date,"_Selected_curve_bandwidth_.tiff"), units="in", width=7, height=7, res=200)
print({
df_selected[,c("labels","LI","LS"):=list(1:100,fitted-1.64*sqrt(0.07189111),
                                         fitted+1.64*sqrt(0.07189111))] %>% 
  ggplot()+
  geom_line(aes(labels,fitted))+geom_line(aes(labels,LI),linetype="twodash")+
  geom_line(aes(labels,LS),linetype="twodash")+ylab("Absorbencia (Curvas seleccionada)")+
  xlab("Medición")+theme_bw()
})
dev.off()



mean_fd_altos=mean.fd(fit_fd_altos$fd)
trimmed_mean_fd_altos=mean.fd(fit_fd_altos$fd,trim=0.1)
std_fd_altos=std.fd(fit_fd_altos$fd)

mean_fd_bajos=mean.fd(fit_fd_bajos$fd)
trimmed_mean_fd_bajos=mean.fd(fit_fd_bajos$fd,trim=0.1)
std_fd_bajos=std.fd(fit_fd_bajos$fd)

summary_stats<-data.frame(
  rbind(
    data.frame(Medicion=1:100,Mean=predict(mean_fd_bajos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Bajo (Grasa<=20)"),
    data.frame(Medicion=1:100,Mean=predict(mean_fd_altos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Alto (Grasa>20)")))

tiff(paste0(date,"_Mean_sd_groups.tiff"), units="in", width=7, height=7, res=200)
print({
ggplot(absorp_4)+
  geom_line(data=absorp_4 ,aes(x=Medicion,y=Fitted,group=Curva),color="darkorange",size=1,alpha=0.14)+
  geom_line(data=summary_stats,aes(x=Medicion,y=mean),size=2,color="darkred",linetype="twodash")+
  geom_line(data=summary_stats,aes(x=Medicion,y=Sd),size=2,color="green",linetype="twodash")+
  guides(color=FALSE)+theme_bw()+xlab("Medición")+ylab("Absorbencia (Curvas ajustadas)")+
  facet_grid(.~Grupo)
})
dev.off()


mean_fd_altos=mean.fd(fit_fd_altos$fd)
trimmed_mean_fd_altos=mean.fd(fit_fd_altos$fd,trim=0.1)
std_fd_altos=std.fd(fit_fd_altos$fd)

mean_fd_bajos=mean.fd(fit_fd_bajos$fd)
trimmed_mean_fd_bajos=mean.fd(fit_fd_bajos$fd,trim=0.1)
std_fd_bajos=std.fd(fit_fd_bajos$fd)

summary_stats_2<-data.frame(
  rbind(
    data.frame(Medicion=1:100,Mean=predict(trimmed_mean_fd_bajos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Bajo (Grasa<=20)"),
    data.frame(Medicion=1:100,Mean=predict(trimmed_mean_fd_altos,newdata=1:100),Sd=predict(std_fd_bajos,newdata=1:100),Grupo="Alto (Grasa>20)")))

tiff(paste0(date,"_Trimmed_mean_sd_per_groups.tiff"), units="in", width=7, height=7, res=200)
print({
ggplot(absorp_4)+
  geom_line(data=absorp_4 ,aes(x=Medicion,y=Fitted,group=Curva),color="darkorange",size=1,alpha=0.14)+
  geom_line(data=summary_stats_2,aes(x=Medicion,y=mean),size=2,color="darkred",linetype="twodash")+
  geom_line(data=summary_stats_2,aes(x=Medicion,y=Sd),size=2,color="green",linetype="twodash")+
  guides(color=FALSE)+theme_bw()+xlab("Medición")+ylab("Absorbencia (Curvas ajustadas)")+
  facet_grid(.~Grupo)
})
dev.off()



#Derivadas

deriv_altos_1<-deriv.fd(fit_fd_altos$fd,1)
deriv_altos_1.1<-predict(deriv_altos_1,1:100)  %>% as.data.table() 
names(deriv_altos_1.1)<-paste0("X",1:ncol(deriv_altos_1.1))



deriv_bajos_1<-deriv.fd(fit_fd_bajos$fd,1)
deriv_bajos_1.1<-predict(deriv_bajos_1,1:100)  %>% as.data.table() 
names(deriv_bajos_1.1)<-paste0("X",1:ncol(deriv_bajos_1.1))


deriv_altos_2<-deriv.fd(fit_fd_altos$fd,2)
deriv_altos_2.1<-predict(deriv_altos_2,1:100)  %>% as.data.table() 
names(deriv_altos_2.1)<-paste0("X",1:ncol(deriv_altos_2.1))


deriv_bajos_2<-deriv.fd(fit_fd_bajos$fd,2)
deriv_bajos_2.1<-predict(deriv_bajos_2,1:100) %>% as.data.table() 
names(deriv_bajos_2.1)<-paste0("X",1:ncol(deriv_bajos_2.1))

derivadas_1<-rbind(deriv_altos_1.1 %>% mutate(labels=1:100) %>% 
                     melt(id.vars="labels",variable.name="Variable",value.name="Valor") %>% mutate(Grupo="Alto (Grasa>20)"),
                   deriv_bajos_1.1 %>% mutate(labels=1:100) %>% 
                     melt(id.vars="labels",variable.name="Variable",value.name="Valor") %>% mutate(Grupo="Bajo (Grasa<=20)"))

derivadas_2<-rbind(deriv_altos_2.1 %>% mutate(labels=1:100) %>% 
                     melt(id.vars="labels",variable.name="Variable",value.name="Valor") %>% mutate(Grupo="Alto (Grasa>20)"),
                   deriv_bajos_2.1 %>% mutate(labels=1:100) %>% 
                     melt(id.vars="labels",variable.name="Variable",value.name="Valor") %>% mutate(Grupo="Bajo (Grasa<=20)"))

tiff(paste0(date,"_Primera_derivada.tiff"), units="in", width=7, height=7, res=200)
print({
  derivadas_1 %>% ggplot()+geom_line(aes(labels,Valor,group=Variable),color="midnightblue")+
    guides(color=F)+theme_bw()+xlab("Medición")+ylab("D(Absorbancia)")+facet_grid(.~Grupo)
})
dev.off()

tiff(paste0(date,"_Segunda_derivada.tiff"), units="in", width=7, height=7, res=200)
print({
derivadas_2 %>% ggplot()+geom_line(aes(labels,Valor,group=Variable),color="midnightblue")+
  guides(color=F)+theme_bw()+xlab("Medición")+ylab("D^2(Absorbancia)")+facet_grid(.~Grupo)
})
dev.off()


W.pca$values
W.pca = pca.fd(fit_selected_curve$fd, nharm=1)
plot(W.pca$harmonics, lwd=3)
sum(W.pca$varprop)
class(W.pca)

tiff(paste0(date,"_Componente_principal.tiff"), units="in", width=7, height=7, res=200)
print({
data.table(labels=1:100,Valor=predict(W.pca$harmonics,newdata=1:100)) %>% 
  ggplot()+geom_line(aes(labels,Valor.PC1),color="midnightblue",size=1.2)+
  ylab("Componente principal")+xlab("Medición")+theme_bw()
})
dev.off()
as.numeric(round(W.pca$harmonics$coefs,2))

