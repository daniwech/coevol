####################################################################################################################
#
library(Matrix)
library(deSolve)

#initial assignments
n=2
k=2
N=c(0.5,0.6)#runif(n)
V=c(0.4,0.5,0.6,0.7)#runif(n*k)
rho=1
lamda=1
eta=1
gamma=0
M=1
r=c(0.5,0.5)#runif(n)
beta=1
c=c(1,1,1,1)#runif(n*k)*10
h=0.5
#A&B together

AB<-matrix(c(
  0,0,0,1,
  1,0,0,1,
  1,0,0,0,
  1,1,1,0),nrow=4,ncol=4,byrow=T
)

if (gamma != 0) {
  cat("WARNING: gamma should be 0")
}

 
#function to Vp-Vq
Vdiff<-function(V){
  k<-length(V)
  out<-matrix(NA,k,k)
  for (p in 1:k){
    for (q in 1:k){
      out[p,q]<-V[p]-V[q]
    }
  }
  return(out)
}

#function to plot
plotout<-function(out,sp_N=2){
  par(mfrow=c(2,2),mar=c(2,2,2,2),mgp=c(1.2,0.2,0),tcl=0.3)
  for(i in 1:(sp_N+1)){
    matplot(out[,1],out[,seq(sp_N*(i-1)+2,sp_N*(i)+1,1)],type="l",ylab="value",xlab="time",main=ifelse(i==1,"population",paste("trait in sp",i-1)))
  }
  legend("topright", legend = 1:sp_N, col = 1:sp_N, lty = 1:sp_N,cex = 0.75)
}

#function of change rate of population N
DN<-function(N,V,AB,n,r,k,beta,gamma,fitness){
  mat <- function(wholesize, blocksize) {
    (matrix(c(rep(1, blocksize), rep(0, wholesize)), wholesize, wholesize/blocksize))
  }
  
  eco_A<-(exp(-(Vdiff(V))^2))*AB
  
  eco_A<-t(mat(k*n,k))%*%eco_A%*%mat(k*n,k)
  
  #diag(eco_A)<--beta
  # Diagonal needs to be zero otherwise, traits within a species
  # also interact with each other (hence have effect on growth rate).
  diag(eco_A)<--0
  DN_out<-N*(r+gamma-beta*N+eco_A%*%(N/(1+h*N)))
  DN_out

}

##function of change rate of traits DV
DV<-function(AB,N,V,c,n,k,M,eta,rho,lamda){
  
  diag_block<-as.matrix(bdiag(replicate(n,list(matrix(1, k,n)))))
  AB_epi<-diag_block*AB
  diag(AB_epi)<--1# get the Vp in to matrix
  mu_self<-c+AB_epi%*%V
  diag(AB_epi)<-0
  epi_contribution<-2*rho*M*(k^-1)*(mu_self/(1+rho*mu_self^2)^2-t(t(mu_self/(1+rho*mu_self^2)^2)%*%AB_epi))
  Vp_q<-Vdiff(V)
  AB_eco<-(diag_block*-1+1)*AB
  eco_contribution<-(-AB_eco*Vp_q*exp(-Vp_q^2))%*%rep(N,each=k)
  DV_out<-lamda*rep(N,each=k)*(epi_contribution+eta*eco_contribution)
  fitness0<-M*(k^-1)*(mu_self/(1+rho*mu_self^2))
  fitness<-as.numeric( tapply(fitness0,cut(1:(n*k),k),FUN=sum))
  list(DV_out=DV_out,fitness=fitness)
}

# model for ODE
model<-function(t,y,parms){
  N=y[1:n]
  V=y[(n+1):(n+(n*k))]
  with(parms,{
    dV=DV(AB=AB,N=N,V=V,c=c,n=n,k=k,rho=rho,lamda=lamda,eta=eta,M=M)
    fitness<-dV$fitness
    dN=DN(N=N,V=V,AB=AB,n=n,r=r,k=k,beta=beta,gamma=gamma,fitness=fitness)
    #y=c(N,V)
    return(list (c(dN,dV$DV_out)))
  })
  
}

#combine parameters
parms = list(AB=AB,r=r,c=c,beta=beta,lamda=lamda,rho=rho,eta=eta,k=k,n=n,M=M,gamma=gamma)

##integration
out <- ode(func = model, y =c(N,V), times = seq(0, 20,0.1), parms =parms)

#plot the results
plotout(out)


#2017.4.21#############################################################################################################
#######################################################################################################################
#######################################################################################################################
library(Matrix)
library(deSolve)

#function to Vp-Vq
Vdiff<-function(V){
  k<-length(V)
  out<-matrix(NA,k,k)
  for (p in 1:k){
    for (q in 1:k){
      out[p,q]<-V[p]-V[q]
    }
  }
  return(out)
}

#function to plot
plotout<-function(out,sp_N=2){
  par(mar=c(2,2,2,2),mgp=c(1.2,0.2,0),tcl=0.3)
  for(i in 1:(sp_N+1)){
    matplot(out[,1],out[,seq(sp_N*(i-1)+2,sp_N*(i)+1,1)],type="l",ylab="value",xlab="time",main=ifelse(i==1,"population",paste("trait in sp",i-1)),ylim=c(0,2))
  }
  legend("topright", legend = 1:sp_N, col = 1:sp_N, lty = 1:sp_N,cex = 0.75)
}

#function of change rate of population N
DN<-function(N,V,AB,n,r,k,beta,gamma,fitness){
  mat <- function(wholesize, blocksize) {
    suppressWarnings(matrix(c(rep(1, blocksize), rep(0, wholesize)), wholesize, wholesize/blocksize))
  }
  
  eco_A<-(exp(-(Vdiff(V))^2))*AB
  eco_A<-t(mat(k*n,k))%*%eco_A%*%mat(k*n,k)
  #diag(eco_A)<--beta
  DN_out<-N*(r+gamma*0-beta*N+eco_A%*%(N/(1+h*N)))
  DN_out
}

##function of change rate of traits DV
DV<-function(AB,N,V,c,n,k,M,eta,rho,lamda,h){
  
  diag_block<-as.matrix(bdiag(replicate(n,list(matrix(1, k,n)))))
  AB_epi<-diag_block*AB
  diag(AB_epi)<--1# get the Vp in to matrix
  mu_self<-c+AB_epi%*%V
  diag(AB_epi)<-0
  epi_contribution<-2*rho*M*(k^-1)*(mu_self/(1+rho*mu_self^2)^2-t(t(mu_self/(1+rho*mu_self^2)^2)%*%AB_epi))
  Vp_q<-Vdiff(V)
  AB_eco<-(diag_block*-1+1)*AB
  eco_contribution<-(-AB_eco*Vp_q*exp(-Vp_q^2))%*%rep(N,each=k)
  DV_out<-lamda*rep(N,each=k)*(epi_contribution+eta*eco_contribution)
  fitness0<-M*(k^-1)*(mu_self/(1+rho*mu_self^2))
  fitness<-as.numeric( tapply(fitness0,cut(1:(n*k),k),FUN=sum))
  list(DV_out=DV_out,fitness=fitness)
}

# model for ODE
model<-function(t,y,parms){

  with(parms,{
    N=y[1:n]
    V=y[(n+1):(n+(n*k))]
    dV=DV(AB=AB,N=N,V=V,c=c,n=n,k=k,rho=rho,lamda=lamda,eta=eta,M=M,h=h)
    fitness<-dV$fitness
    dN=DN(N=N,V=V,AB=AB,n=n,r=r,k=k,beta=beta,gamma=gamma,fitness=fitness)
    #y=c(N,V)
    return(list (c(dN,dV$DV_out)))
  })
  
}

#function to summary results of each simulation
eco_evo22<-function(m,parms,N,V){
  #parms<-c(AB=list(m),parms0)
  parms$AB<- m
  out <-ode(func = model, y =c(N,V), times = seq(0, 1000), parms =parms)
  
  #abundance at the end of simulation
  abun<-as.numeric(t(tail(out[,-1],1)))
  #variance
  variance<-apply(tail(out[,-1],100),2,var)<=1e-5
  #the direction of last 100 steps
  slope<-apply(tail(out[,-1],100),2,function(x){round(lm(x~c(1:100))$coefficients[2],2)})
  #time to reach equilibrium;smallest steps when abundance==equilibrium on the condition of variance is 0
  time2eq<-sapply(1:length(abun),function(x){ifelse(variance[x],min(which(abs(out[,1+x]-abun[x])<=1e-5,T)),NA)})
  #wrap them
  return(data.frame(abun=abun,variance,slope,time2eq))
  #plot(0,frame.plot = F,axes=F,type="l",ylab="",xlab="",ylim=c(0,1),xlim=c(0,1))
  #plotrix::addtable2plot (x=0.5,y=0,table=AB_l[[2]],ypad=0.1)
  #plotout(out)
}

#initial assignments
n=2
k=2
N=c(0.5,0.5)#runif(n)
V=c(0.5,0.5,0.5,0.5)#runif(n*k)
rho=1
lamda=0.1
eta=1
gamma=0
M=1
r=c(0.5,0.5)#runif(n)
beta=1
c=c(1,1,1,1)#runif(n*k)*10
h=0.5
parms = list(r=r,c=c,beta=beta,lamda=lamda,rho=rho,eta=eta,k=k,n=n,M=M,gamma=gamma,h=h)

# construction different types of B matrix 
Bs<-matrix(c(0,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,1,1,0,
             0,-1,0,0,
             0,0,-1,0,
             0,-1,-1,0,
             0,-1,1,0,
             0,1,-1,0),ncol=4,byrow=T)
Bs_l<-apply(Bs,1,function(x){list(matrix(x,2,2,byrow=T))})
Bs_l<-unlist(Bs_l,F)

# construction different types of A matrix 
As<-as.matrix(expand.grid(c(0,1,-1),c(0,1,-1),c(0,1,-1),c(0,1,-1)))
As_l<-apply(As,1,function(x){list(matrix(x,2,2,byrow=T))})
As_l<-unlist(As_l,F)

# construction different types of AB matrix by combining A and B matrix
AB_l<-apply(expand.grid(c(1:9),c(1:9),c(1:81),c(1:81)),1,function(x){list(rbind(cbind(Bs_l[[x[1]]],As_l[[x[3]]]),cbind(As_l[[x[4]]],Bs_l[[x[2]]])))})
AB_l<-unlist(AB_l,F)

#run the simulation parallelly
Sys.time()
library(snowfall)
sfInit(parallel=TRUE, cpus=46)
sfLibrary(deSolve)
sfLibrary(Matrix)
sfExport('Vdiff','DN','DV','model','eco_evo22','AB_l','parms','N','V','h')

res<-sfLapply(1:length(AB_l),function(i){tryCatch(eco_evo22(AB_l[[i]],parms=parms,N=N,V=V),error=function(e) NULL)})
sfStop()
Sys.time()

save.image(file = "D:/work/eco_evo/2017.2.11.RData")

##############################################################
### summary of the structure of AB matrix on the coexistence##
### of population and traits and both                       ##
##############################################################
#B1 matrix to locate inner and inter part of AB
 B1<-matrix(c(1,1,0,0,
         1,1,0,0,
         0,0,1,1,
         0,0,1,1),nrow=4,byrow=T)==1

lower.loc<-which(lower.tri(B1)==T,arr.ind=T)#lowerpart 
upper.loc<-lower.loc[,c(2,1)]# and upper 


library(ggplot2);library(cowplot)
#sturcture of AB: total_C, connectance of AB; B_C, connectance of B; A_C, connectance of A;A_Positive, sum of AB; B_positive, sum of AB; symm, sum of (lower-upper)
M_structure<-t(sapply(AB_l,function(x){c(total_C=sum(x!=0),B_C=sum(x[B1]!=0),A_C=sum(x[!B1]!=0),A_Positive=sum(x[!B1]),B_positive=sum(x[B1]), symm=sum(x[lower.loc]-x[upper.loc]))}))

#coexistence of population (abun_coexist), trait (trait_coexist) and both (all_coexist)
coexist<-t(sapply(res,function(x){c(abun_coexist=all(x[1:2,1]>1e-5),trait_coexist=all(x[3:6,1]>1e-5),all_coexist=all(x[1:6,1]>1e-5))}))
coexist<-apply(coexist,2,as.numeric)

#gam spline function to show the trend of coexistence against structure of AB
tmp<-data.frame(coexist,M_structure)
  #tmp$coexist<-as.factor(tmp$coexist)     
ps.1<-lapply(colnames(M_structure),function(x){ggplot(tmp,aes_string(x,"abun_coexist")) + geom_smooth( method="gam",formula=y~s(x,bs="cr",k=4),method.args = list(family = "binomial"))})
ps.2<-lapply(colnames(M_structure),function(x){ggplot(tmp,aes_string(x,"trait_coexist")) + geom_smooth( method="gam",formula=y~s(x,bs="cr",k=4),method.args = list(family = "binomial"))})
ps.3<-lapply(colnames(M_structure),function(x){ggplot(tmp,aes_string(x,"all_coexist")) + geom_smooth( method="gam",formula=y~s(x,bs="cr",k=4),method.args = list(family = "binomial"))})

#save plots
ggsave("D:/work/eco_evo/plots/ps.1abun_coexist.png",cowplot::plot_grid(plotlist=ps.1, ncol=2,nrow=3,align="hv"),width=16, height=16,units="cm")
ggsave("D:/work/eco_evo/plots/ps.2trait_coexist.png",cowplot::plot_grid(plotlist=ps.2, ncol=2,nrow=3,align="hv"),width=16, height=16,units="cm")
ggsave("D:/work/eco_evo/plots/ps.3all_coexist.png",cowplot::plot_grid(plotlist=ps.3, ncol=2,nrow=3,align="hv"),width=16, height=16,units="cm")


