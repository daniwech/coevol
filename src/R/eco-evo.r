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
  0,1,1,1,
  1,0,1,1,
  1,1,0,1,
  1,1,1,0),nrow=4,ncol=4,byrow=T
)

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
    suppressWarnings(matrix(c(rep(1, blocksize), rep(0, wholesize)), wholesize, wholesize/blocksize))
  }
  
  eco_A<-(exp(-(Vdiff(V))^2))*AB
  eco_A<-t(mat(k*n,k))%*%eco_A%*%mat(k*n,k)
  #diag(eco_A)<--beta
  DN_out<-N*(r+gamma*0-beta*N+eco_A%*%(N/(1+h*N)))
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
