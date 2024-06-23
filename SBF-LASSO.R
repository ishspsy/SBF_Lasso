library(SAM)
library(foreach)
library(caTools)
integ<-trapz
library(tmvtnorm)
library(truncnorm)
library(doSNOW)
library(doMC)
library(tictoc)


###
### construction of the normalized kernel
###

# Epanechnikov kernel
K<-function(t)  (1-t^2)*3/4*dunif(t,-1,1)*2


# scaling of the kernel
Kh<-function(u,v,h)    K((u-v)/h)/h

# normalized kernel
K_h<-function(s,t,g){
  numer<-Kh(s,t,g)
  
  if(t<a || t>b){
    return(rep(0,length(s)))
  }else{
    denom<-integrate(Kh,lower=a,upper=b,v=t,h=g)$val
    
    if(denom==0){
      return(rep(0,length(s)))
    }else{
      return(numer/denom)
    }
  }
}



###
### kernel density estimation (joint and marginal)
###


p_0<-function(X){
  tmp<-rep(1,nrow(X))
  for(j in 1:ncol(X)){
    tmp<-tmp*dunif(X[,j],a,b)*(b-a)}
  return(mean(tmp))}

p_kj<-function(x_k, x_j, X, ind_kj, h, K_h){
  # x_k: evaluation points of the k th component
  # x_j: evaluation points of the j th component
  # X: n0 by d0 covariate matrix
  # ind: 2-dimensional index vector (e.g.) ind=c(k,j)
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  
  n0<-nrow(X)
  k<-ind_kj[1]
  j<-ind_kj[2]
  
  p_hat<-c()
 
  
  registerDoMC(cores=3)
  tmp=0
  tmp=foreach(i = 1:n0, .combine='+')  %dopar% NW0(i,j,k,a,b,x_k,X,h,K_h,x_j)   

  p_hat<-tmp/n0
  return(p_hat)}


NW0 <- function(i,j,k,a,b,x_k,X,h,K_h,x_j){
  K_h(x_k,X[i,k],h[k])%*%t(K_h(x_j,X[i,j],h[j]))*prod(dunif(X[i,],a,b))*(b-a)^(length(X[i,]))}


NW1 <- function(i,j,x,X,h,Y,a,b,K_h){
  K_h(x[,j],X[i,j],h[j])*Y[i]*prod(dunif(X[i,],a,b))*(b-a)^(length(X[i,]))
}

NW2 <- function(i,j,x,X,h,Y,a,b,K_h){
  K_h(x[,j],X[i,j],h[j])*prod(dunif(X[i,],a,b))*(b-a)^(length(X[i,]))
}

NW3 <- function(i,j,x,X,h,Y,a,b,x_j,K_h){
  K_h(x_j,X[i,j],h[j])*prod(dunif(X[i,],a,b))*(b-a)^(length(X[i,]))
}


p_j<-function(x_j, X, ind_j, h, K_h){
  # x_j: evaluation points of the j th component
  # X: n0 by d0 covariate matrix
  # ind: index scalar (e.g.) ind=j
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  
  n0<-nrow(X)
  j<-ind_j
  
  registerDoMC(cores=3)
  tmp=0
  tmp=foreach(i = 1:n0, .combine='+')  %dopar% NW3(i,j,x,X,h,Y,a,b,x_j,K_h)    
  
  p_hat<-tmp/n0
  return(p_hat)}



###
### projection operator of the k th component function to the j th component space
###
proj_kj<-function(x_j, f, X, ind_kj, h, K_h){
  # x_j: evaluation points of the j th component
  # f: matrix for the evaluated values at points for each component function
  # X: n0 by d0 covariate matrix
  # ind: 2-dimensional index vector (e.g.) ind=c(k,j)
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  
  k<-ind_kj[1]
  j<-ind_kj[2]
  
  f_k<-f[,k]
  x_k<-c()
  if(length(f_k)==nrow(X)){
    x_k<-X[,k]
  }else{
    x_k<-x_j}
  
  asdf<-p_j(x_j,X,j,h,K_h)
  tmp_ind<-which(asdf!=0)
  
  if(length(tmp_ind)>0){
    p_hat<-matrix(0,nrow=length(x_k),ncol=length(x_j))
    p_hat[,tmp_ind]<-t(t(p_kj(x_k,x_j[tmp_ind],X,c(k,j),h,K_h))/asdf[tmp_ind])
    
    tmp<-c()
    for(l in 1:ncol(p_hat)){
      tmptmp<-f_k*c(p_hat[,l])
      tmp[l]<-integ(sort(x_k),tmptmp[order(x_k)])
    }
    
    return(tmp)    
    
  }else{
    return(0)
  }
}



###
### Nadaraya-Watson estimator
###
NW<-function(Y,X,x,h,K_h){
  # Y: n0-dimensional response vector
  # X: n0 by d0 covariate matrix
  # x: N by d0 covariate matrix for evaluation
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  
  n0<-length(Y)
  d0<-ncol(X)
  
  N<-nrow(x)
  
  f_nw<-matrix(0,nrow=N,ncol=d0)
  for(j in 1:d0){
    registerDoMC(cores=3)
    tmp1=foreach(i = 1:n0, .combine='+') %dopar% NW1(i=i,j=j,x=x,X=X,h=h,Y=Y,a=a,b=b,K_h=K_h)   
    registerDoMC(cores=3)
    tmp2=foreach(i = 1:n0, .combine='+') %dopar% NW2(i=i,j=j,x=x,X=X,h=h,Y=Y,a=a,b=b,K_h=K_h) 
    
    tmp_ind<-which(tmp2!=0)
    f_nw[tmp_ind,j]<-tmp1[tmp_ind]/tmp2[tmp_ind]}
  return(f_nw)
}


###
### smooth backfitting algorithm
###

cent<-function(x,X,j,h,K_h,f){
  tmp1<-p_j(x,X,j,h,K_h)
  tmp<-f-integ(sort(x),(f*tmp1)[order(x)])
  return(tmp)
}


SBF_update_j<-function(Y,X,x,f,f0,ind_up,h,K_h,f_nw){
  # Y: n0-dimensional response vector
  # X: n0 by d0 covariate matrix (obs)
  # x: N by d0 covarate matrix (evaluation point)
  # f: N by d0 component function matrix (evaluated at x) to be updated
  # f0: N by d0 component function matrix (evaluated at x) to be used in updating f
  # ind_up: (scalar) index of component to be updated
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  # f_nw: N by d0 matrix of Nadaraya-Watson estimator for component functions
  
  tmp_index<-which(apply(apply(X,1,'dunif',min=a,max=b),2,'prod')>0)
  Y_mean<-sum(Y[tmp_index])/length(Y)/p_0(X)
  
  j<-ind_up
  
  tmp1<-tmp2<-0
  
  if(j==1){
    registerDoMC(cores=3)
    tmp2=foreach(k = (j+1):ncol(X), .combine='+')  %dopar% proj_kj(x[,j],f0,X,c(k,j),h,K_h)
    f[,j]<-f_nw[,j]-Y_mean-tmp2 
  }
  if(j>1 && j<ncol(X)){
    
    registerDoMC(cores=3)
    tmp1=foreach(k = 1:(j-1), .combine='+')  %dopar% proj_kj(x[,j],f,X,c(k,j),h,K_h)

    registerDoMC(cores=3)
    tmp2=foreach(k = (j+1):ncol(X), .combine='+')  %dopar% proj_kj(x[,j],f0,X,c(k,j),h,K_h)

    f[,j]<-f_nw[,j]-Y_mean-tmp1-tmp2
  }
  if(j==ncol(X)){
    registerDoMC(cores=3)
    tmp2=foreach(k = 1:(ncol(X)-1), .combine='+')  %dopar% proj_kj(x[,ncol(X)],f,X,c(k,ncol(X)),h,K_h)
    f[,ncol(X)]<-f_nw[,ncol(X)]-Y_mean-tmp2 
  }
  f[,j]<-cent(x[,j],X,j,h,K_h,f[,j])
  return(f)
}


shrink_update_j<-function(Y,X,x,f,ind_up,h,K_h,lambda){
  # Y: n0-dimensional response vector
  # X: n0 by d0 covariate matrix (obs)
  # x: N by d0 covarate matrix (evaluation point)
  # f: N by d0 component function matrix (evaluated at x) to be shrinken
  # ind_up: (scalar) index of component to be updated
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  # f_nw: N by d0 matrix of Nadaraya-Watson estimator for component functions
  # lambda: (scalar) thresholding parameter
  
  j<-ind_up
  
  x_j<-c(x[,j])
  
  p_hat_j<-p_j(x_j,X,j,h,K_h)
  sq_f_j<-c(f[,j]^2)*c(p_hat_j)
  
  norm_f_j<-sqrt(integ(sort(x_j),sq_f_j[order(x_j)]))
  
  
  if(norm_f_j==0){
    f[,j]<-0
  }else{
    
    if(norm_f_j>lambda){
      f[,j]<-(1-lambda/norm_f_j)*f[,j]
    }else{
      f[,j]<-0
    }
  }
  f[,j]<-cent(x[,j],X,j,h,K_h,f[,j])
  #f[,j] <- f[,j] - mean(f[,j])
  return(f)  
}



# SBF algorithm
sparse_SBF_fit<-function(Y,X,x,h,K_h,f_nw,lambda,fini){
  # Y: n0-dimensional response vector
  # X: n0 by d0 covariate matrix (obs)
  # x: N by d0 covarate matrix (evaluation point)
  # h: d0-dimensional bandwidth vector
  # K_h: (univariate) function class object corresponding to the normalized kernel function
  # f_nw: N by d0 matrix of Nadaraya-Watson estimator for component functions
  # lambda: (scalar) thresholding parameter
  
  n0<-length(Y)
  d0<-ncol(X)
  
  N<-nrow(x)
  
  tmp_index<-which(apply(apply(X,1,'dunif',min=a,max=b),2,'prod')>0)
  Y_mean<-sum(Y[tmp_index])/length(Y)/p_0(X)
  
  # initial estimator
  f<-matrix(nrow=N,ncol=d0)
  for(j in 1:d0){
    f[,j]<-cent(x[,j],X,j,h,K_h,f_nw[,j])
  }
  
  f = fini
  
  # backfitting
  eps<-10
  tmp_eps<-20
  iter<-1
  while(eps>0.2){
    #print(eps)
    tmp_eps<-eps
    f0<-f
    f1<-f
    ff1<-f

        for(j in 1:d0){
          f[,j]<-SBF_update_j(Y,X,x,f,f,j,h,K_h,f_nw)[,j]
          f[,j]<-shrink_update_j(Y,X,x,f,j,h,K_h,lambda)[,j]}

    eps<-max(sqrt(apply(abs(f-f0)^2,2,'mean')))
    

    if(iter>30){
      return(f0)
    }
    iter<-iter+1
  }
  return(f)
}





### MC simulation

sim<-function(M,n0,d0,h,lambda,corr,esd, case, fini){
  
  f_fit<-matrix(0,nrow=N,ncol=d0)
  zero_call<-rep(0,d0)
  
  m=M
  set.seed(m)
  print(m)
  
  if (corr>0.01){
    S = array(0,c(d0,d0))
    for (i1 in 1:d0){for (j1 in 1:d0){
      S[i1,j1] = corr^(abs(i1-j1))}}
    
    X<-rtmvnorm(n0,sigma=S,lower=rep(a,d0),upper=rep(b,d0))}
  
  if (case==1 | case==3){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)}
  
  if (case==2 | case==4){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)
    U= matrix(runif(n0,0,1), nrow=n0, ncol=1)
    U = matrix(rep(U,d0), nrow=n0, byrow=F)
    X= 0.5*X+0.5*U}
  
  
  f_tr<-matrix(0,nrow=n0,ncol=d0)
  for(i in 1:n0){
    if (case <= 2){
      f_tr[i,1]<-f1(X[i,1])
      f_tr[i,2]<-f2(X[i,2])
      f_tr[i,3]<-f3(X[i,3])
      f_tr[i,4]<-f4(X[i,4])}
    if (case >= 3){
      f_tr[i,1]<-f1(X[i,1])
      f_tr[i,2]<-f2(X[i,2])
      f_tr[i,3]<-f3(X[i,3])
      f_tr[i,4]<-f4(X[i,4])
      f_tr[i,5]<-f5(X[i,5])
      f_tr[i,6]<-f6(X[i,6])
      f_tr[i,7]<-f7(X[i,7])
      f_tr[i,8]<-f8(X[i,8])
      f_tr[i,9]<-f9(X[i,9])
      f_tr[i,10]<-f10(X[i,10])
      f_tr[i,11]<-f11(X[i,11])
      f_tr[i,12]<-f12(X[i,12])
    }}
  
  Y<-c()
  for(i in 1:n0){
    Y[i]<-sum(f_tr[i,])+rnorm(1,0,esd)
  }
  
  f_nw<-NW(Y,X,x,h,K_h) 
  tic()
  f_tmp<-sparse_SBF_fit(Y,X,x,h,K_h,f_nw,lambda,fini)
  toc1=toc()
  norm_tmp<-apply(f_tmp^2*diff(x0)[1],2,'sum')
  print(round(norm_tmp,3))
  zero_call<-zero_call+(norm_tmp<0.005)
  
  f_fit<- f_tmp
  
  return(list(zero=zero_call,f=f_fit,toc1=toc1))
}


interpol = function(x, f_fit, a, b){
  p = length(x)
  points_ab = seq(a, b, length.out=101)
  values = 0
  for (j in 1:p){
    xj = x[j]; 
    if (abs(xj)<0.000001){xj=0.000001}
    f_fitj = f_fit[,j]
    cutxj = as.numeric(cut(xj, breaks=points_ab, labels= 1:100))
    a1 = points_ab[cutxj]; b1 = points_ab[cutxj+1]
    ya1 = f_fitj[cutxj]; yb1 = f_fitj[cutxj+1]
    w1 = (xj-a1) / (b1-a1); w2=1-w1
    values = values + ya1 * w2 + yb1 * w1
  }
  return(values)}


h_RT <- function(X, Y) {
  
  # Quartic fit
  mod_Q <- lm(Y ~ poly(X, raw = TRUE, degree = 4))
  
  # Estimates of unknown quantities
  int_sigma2_hat <- diff(range(X)) * sum(mod_Q$residuals^2) / mod_Q$df.residual
  theta_22_hat <- mean((2 * mod_Q$coefficients[3] +
                          6 * mod_Q$coefficients[4] * X +
                          12 * mod_Q$coefficients[5] * X^2)^2)
  # h_RT
  #R_K <- 0.5 / sqrt(pi)
  R_K <- 3/5
  ((R_K * int_sigma2_hat) / (theta_22_hat * length(X)))^(1 / 5)
}




h_opt_ftn<-function(M,n0,d0,corr,esd, case, a, b){
  points_ab = seq(a,b, length.out=N)
  
  zero_call<-rep(0,d0)
  
  m=M
  set.seed(m)
  
  if (corr>0.01){
    S = array(0,c(d0,d0))
    for (i1 in 1:d0){for (j1 in 1:d0){
      S[i1,j1] = corr^(abs(i1-j1))}}
    
    X<-rtmvnorm(n0,sigma=S,lower=rep(a,d0),upper=rep(b,d0))}
  
  if (case==1 | case==3){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)}
  
  if (case==2 | case==4){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)
    U= matrix(runif(n0,0,1), nrow=n0, ncol=1)
    U = matrix(rep(U,d0), nrow=n0, byrow=F)
    X= 0.5*X+0.5*U}
  
  
  f_tr<-matrix(0,nrow=n0,ncol=d0)
  for(i in 1:n0){
    if (case <= 2){
      f_tr[i,1]<-f1(X[i,1])
      f_tr[i,2]<-f2(X[i,2])
      f_tr[i,3]<-f3(X[i,3])
      f_tr[i,4]<-f4(X[i,4])}
    if (case >= 3){
      f_tr[i,1]<-f1(X[i,1])
      f_tr[i,2]<-f2(X[i,2])
      f_tr[i,3]<-f3(X[i,3])
      f_tr[i,4]<-f4(X[i,4])
      f_tr[i,5]<-f5(X[i,5])
      f_tr[i,6]<-f6(X[i,6])
      f_tr[i,7]<-f7(X[i,7])
      f_tr[i,8]<-f8(X[i,8])
      f_tr[i,9]<-f9(X[i,9])
      f_tr[i,10]<-f10(X[i,10])
      f_tr[i,11]<-f11(X[i,11])
      f_tr[i,12]<-f12(X[i,12])
    }}
  
  Y<-c()
  for(i in 1:n0){
    Y[i]<-sum(f_tr[i,])+rnorm(1,0,esd)
  }
  
  h_set = NULL
  for (ii in 1:ncol(X)){
    data_ind = cbind(Y,X[,ii])
    h_set=c(h_set,h_RT(X[,ii],Y))
  }
  
  
  return(h_set)
}



# number of components
d0 = 10
# sample size
n0= 100


# number of evaluation points
N<-101

## Simulation setting
case =1
corr = 0   # covariate correlation
esd=1      # error variance

###
### simulation
###

# covariate range 
a<- 0
b<- 1
t<-seq(a,b,length.out=N)


if (case <= 2){
  # Example 1 and 2
  ### component functions 
  g01<-function(t) -sin(2*t)
  g02<-function(t) t^2-25/12
  g03<-function(t) t
  g04<-function(t) exp(-t)-2/5*sinh(5/2)
  
  
  # Example 1 and 2
  ### component functions
  g01<-function(t) -sin(2*t)
  g02<-function(t) t^2-25/12
  g03<-function(t) t
  g04<-function(t) exp(-t)-2/5*sinh(5/2)
  
  
  # Example 1 and 2
  g01<-function(t) 5*t - 2.5
  g02<-function(t) 3*(2*t-1)^2 -1
  g03<-function(t) 4*sin(2*pi*t) / (2-sin(2*pi*t)) -0.619
  g04<-function(t) 6*(0.1*sin(2*pi*t) + 0.2*cos(2*pi*t)+0.3*sin(2*pi*t)^2+
                        0.4*cos(2*pi*t)^3 + 0.5*sin(2*pi*t)^3) - 0.9
  
  c(integrate(g01, 0,1),integrate(g02, 0,1),integrate(g03, 0,1),integrate(g04, 0,1))
  
  g1<-function(t)    g01(t)*dunif(t,a,b)
  g2<-function(t)    g02(t)*dunif(t,a,b)
  g3<-function(t)    g03(t)*dunif(t,a,b)
  g4<-function(t)    g04(t)*dunif(t,a,b)
  
  f1<-function(t) g01(t)-integrate(g1,lower=a,upper=b)$val
  f2<-function(t) g02(t)-integrate(g2,lower=a,upper=b)$val
  f3<-function(t) g03(t)-integrate(g3,lower=a,upper=b)$val
  f4<-function(t) g04(t)-integrate(g4,lower=a,upper=b)$val}


if(case >= 3){
  
  # Example 4
  g01<-function(t) 1*t - 0.5
  g02<-function(t) 1*(2*t-1)^2 - 0.333
  g03<-function(t) 1*sin(2*pi*t) / (2-sin(2*pi*t)) - 0.155
  g04<-function(t) 1*(0.1*sin(2*pi*t) + 0.2*cos(2*pi*t)+0.3*sin(2*pi*t)^2+
                        0.4*cos(2*pi*t)^3 + 0.5*sin(2*pi*t)^3) - 0.150
  g05<-function(t) 1.5*1 * t - 0.5 * 1.5
  g06<-function(t) 1.5*1*(2*t-1)^2- 0.333 *1.5
  g07<-function(t) 1.5*1*sin(2*pi*t) / (2-sin(2*pi*t))- 0.155*1.5
  g08<-function(t) 1.5*1*(0.1*sin(2*pi*t) + 0.2*cos(2*pi*t)+0.3*sin(2*pi*t)^2+
                            0.4*cos(2*pi*t)^3 + 0.5*sin(2*pi*t)^3)- 0.150*1.5
  g09<-function(t) 2*1*t -  0.5*2
  g010<-function(t) 2*1*(2*t-1)^2 - 0.333*2
  g011<-function(t) 2*1*sin(2*pi*t) / (2-sin(2*pi*t)) - 0.155*2
  g012<-function(t) 2*1*(0.1*sin(2*pi*t) + 0.2*cos(2*pi*t)+0.3*sin(2*pi*t)^2+
                           0.4*cos(2*pi*t)^3 + 0.5*sin(2*pi*t)^3) - 0.150*2
  
  
  f1<-function(t) g01(t)-integrate(g01,lower=a,upper=b)$val
  f2<-function(t) g02(t)-integrate(g02,lower=a,upper=b)$val
  f3<-function(t) g03(t)-integrate(g03,lower=a,upper=b)$val
  f4<-function(t) g04(t)-integrate(g04,lower=a,upper=b)$val
  f5<-function(t) g05(t)-integrate(g05,lower=a,upper=b)$val
  f6<-function(t) g06(t)-integrate(g06,lower=a,upper=b)$val
  f7<-function(t) g07(t)-integrate(g07,lower=a,upper=b)$val
  f8<-function(t) g08(t)-integrate(g08,lower=a,upper=b)$val
  f9<-function(t) g09(t)-integrate(g09,lower=a,upper=b)$val
  f10<-function(t) g010(t)-integrate(g010,lower=a,upper=b)$val
  f11<-function(t) g011(t)-integrate(g011,lower=a,upper=b)$val
  f12<-function(t) g012(t)-integrate(g012,lower=a,upper=b)$val
}


# evaluation points
x0<-seq(a,b,length.out=N)
x<-matrix(rep(x0,d0),nrow=N,ncol=d0)

f_eval<-matrix(0,nrow=N,ncol=d0)
if (case <= 2){
  f_eval[,1]<-f1(x0)
  f_eval[,2]<-f2(x0)
  f_eval[,3]<-f3(x0)
  f_eval[,4]<-f4(x0)}

if (case >= 3){
  f_eval[,1]<-f1(x0)
  f_eval[,2]<-f2(x0)
  f_eval[,3]<-f3(x0)
  f_eval[,4]<-f4(x0)
  f_eval[,5]<-f5(x0)
  f_eval[,6]<-f6(x0)
  f_eval[,7]<-f7(x0)
  f_eval[,8]<-f8(x0)
  f_eval[,9]<-f9(x0)
  f_eval[,10]<-f10(x0)
  f_eval[,11]<-f11(x0)
  f_eval[,12]<-f12(x0)}



par(mfrow=c(2,2))
plot(t,f1(t),type='l',xlab=paste('x',1,sep=''),ylab=paste('f',1,set=''),
     xaxs='i',xlim=c(a,b),lwd=2)
abline(h=0,col=8)
plot(t,f2(t),type='l',xlab=paste('x',2,sep=''),ylab=paste('f',2,set=''),
     xaxs='i',xlim=c(a,b),lwd=2)
abline(h=0,col=8)
plot(t,f3(t),type='l',xlab=paste('x',3,sep=''),ylab=paste('f',3,set=''),
     xaxs='i',xlim=c(a,b),lwd=2)
abline(h=0,col=8)
plot(t,f4(t),type='l',xlab=paste('x',4,sep=''),ylab=paste('f',4,set=''),
     xaxs='i',xlim=c(a,b),lwd=2)
abline(h=0,col=8)


## M: seed number
M=1
fini = array(0, c(N,d0))
h = h_opt_ftn(M=M,n0,d0,corr=corr,esd=esd, case=case, a, b)
result = sim(lambda =0.5, M=M, n0=n0,d0=d0,h=h,corr=corr,esd=esd,case=case, fini=fini)

## Draw coefficient functions
par(mfrow=c(4,3))
for (jj in 1:d0){
  plot(x0,result$f[,jj],col="blue")
  lines(x0,f_eval[,jj],col="red")
}



## Debiased inference

deb_fun = function(M){
  
  ##### here new (1st)
  registerDoMC(cores=4)
  h_tmp1= foreach(i = M, .combine = 'list', .multicombine = T) %dopar% h_opt_ftn(M=i,n0,d0,corr=corr,esd=esd, case=case, a, b)
  h =h_tmp1
  fini = array(0, c(N,d0))
  res_new = sim(lambda = 0.1, M=M,n0=n0,d0=d0,h=h,corr=corr,esd=esd,SCAD=0,case=case, fini=fini,pb=0)
  est_f = res_new$f
  write.csv(est_f, file= paste("est_f_", M, ".csv", sep=""))
  
  
  ##### (2nd) copy the code (inside the sim function) in the note to create X, Y
  f_fit<-matrix(0,nrow=N,ncol=d0)
  zero_call<-rep(0,d0)
  
  m=M
  set.seed(m)
  print(m)
  
  X<-matrix(runif(n0*d0,a,b),nrow=n0,ncol=d0)
  
  if (corr>0.01){
    S = array(0,c(d0,d0))
    for (i1 in 1:d0){for (j1 in 1:d0){
      S[i1,j1] = corr^(abs(i1-j1))}}
    
    X<-rtmvnorm(n0,sigma=S,lower=rep(a,d0),upper=rep(b,d0))}
  
  if (case==3){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)}
  
  if (case==3.5){
    X<-matrix(runif(n0*d0,0,1),nrow=n0,ncol=d0)
    U= matrix(runif(n0,0,1), nrow=n0, ncol=1)
    U = matrix(rep(U,d0), nrow=n0, byrow=F)
    X= (X+U)/2}
  
  f_tr<-matrix(0,nrow=n0,ncol=d0)
  for(i in 1:n0){
    f_tr[i,1]<-f1(X[i,1])
    f_tr[i,2]<-f2(X[i,2])
    f_tr[i,3]<-f3(X[i,3])
    f_tr[i,4]<-f4(X[i,4])
  }
  Y<-c()
  for(i in 1:n0){
    Y[i]<-sum(f_tr[i,])+rnorm(1,0,esd)
  }
  f_nw<-NW(Y,X,x,h,K_h) 
  
  
  # (3rd) prepare optimization function
  ########### here ######
  m=length(x0)
  matpb = matrix(rep(1:d0, d0), ncol=d0); matpb2 =t(matpb)
  parpjkj = function(mm){
    j = matpb2[mm]; i = matpb[mm];
    if (i != j){
      return(pjkj(x0,  X, c(i,j), h, K_h))}
    if (i == j){
      return(array(1, c(m,m)))
    }
  }
  
  registerDoMC(cores=3)
  tmp2=foreach(k = 1:(d0^2), .combine='cbind')  %dopar% parpjkj(k)
  
  phi=array(0, c(m*d0, m*d0))
  jjj=0
  for (mm in 1:(d0^2)){
    jjj=jjj+1
    j = matpb2[mm]; i = matpb[mm];
    indj = (m*(j-1)+1):(m*j);     indi = (m*(i-1)+1):(m*i);
    phi[indj,indi] = tmp2[, (m*(jjj-1)+1) : (m*jjj)]
  }
  
  phiphi = phi; hh=h; est_f_est_f=est_f;  XX=X # be careful when running
  
  
  ####### inference j_true
  
  
  
  #finaljt = function(j_true){ 
  
  phi=phiphi; h=hh; est_f= est_f_est_f; X=XX
  ind1=1:m; ind2= (m*(j_true-1)+1) : (m*j_true)
  indtotal = 1:(m*d0); indtotal[ind1]=ind2; indtotal[ind2]=ind1; 
  phi = phiphi[indtotal,indtotal]
  h[j_true]=hh[1]; h[1]=hh[j_true];
  est_f[,j_true] = est_f_est_f[,1]; est_f[,1] = est_f_est_f[,j_true]
  X[,j_true]= XX[,1]; X[,1]=XX[,j_true]
  
  write.csv(phi, file=paste("phidata_", M, "_", j_true, ".csv", sep=""))
  
  matlab.lines <- c(
    "d= 10;",
    "m = 101;",
    "ith=[ones(m,m), zeros(m, m*(d-1))];",
    "ith2 = zeros(m*d,d);",
    "for i=1:d",
    "indd = (m*(i-1)+1) : (m*i);",
    "ith2(indd,i)=1;",
    "end",
    "phi = zeros(m*d,m*d);",
    "for i=1:d",
    "indd = (m*(i-1)+1) : (m*i);",
    "phi(indd,indd) = 1+phi(indd,indd);",
    "end",
    "gamma=0.02;",
    "aaa=mfilename('fullpath')",
    "aaaind=strfind(aaa, '_');",
    "aaaind_length= length(aaaind); aaaind1 = aaaind(aaaind_length-1); aaaind2 = aaaind(aaaind_length);",
    "M_ind = aaa((aaaind1+1):(aaaind2-1)); M = str2double(M_ind);",
    "j_ind = aaa((aaaind2+1):(length(aaa))); j_true = str2double(j_ind);",
    "phi=csvread(sprintf('phidata_%d_%d.csv', M, j_true),1,1);",
    "cvx_begin", 
    "variable Q(m, m, d)",
    "minimize(sum(max(norms(Q,2,2))))",
    "subject to",
    "Q(:,:,1)== 0",
    "((ith-(Q(:,:)+ ith)*phi/m).^2) *  ith2  <= m*gamma^2",
    "cvx_end",
    "csvwrite(sprintf('Q_%d_%d.csv', M, j_true),Q(:,:),1,0)"
  )
  date2 = paste("opt_smooth_", M, "_", j_true, ".m", sep="")
  writeLines(matlab.lines, con=date2)
  
  ith=cbind(array(1,c(m,m)), array(0,c(m, m*(d0-1))));
  ithh = array(0, c(m*d0,m*d0))
  
  ith2 = array(0,c(m*d0,d0));
  for (i in 1:d0){
    indd = (m*(i-1)+1) : (m*i);
    ith2[indd,i]=1;
  }
  
  
  # (4th) after optimization function
  Qest = read.csv(paste("Q_", M, "_", j_true, ".csv", sep=""), header=TRUE)   # phi:  dm by dm
  Qest=as.matrix(Qest)
  
  est_f_debiased = matrix(est_f[,1], ncol=1)+ (Qest +ith) %*% (matrix(f_nw,ncol=1) - matrix(rep(mean(Y),m*d0),ncol=1) - phi %*% matrix(est_f,ncol=1)/m)/m^2
  
  serror_fun =function(i){
    Khh = NULL
    for (jj in 1:d0){
      Khh = rbind(Khh, matrix(K_h(x0,X[i,jj],h[jj]),ncol=1))
    }
    
    pjx=p_j(x0, X, 1, h, K_h)
    return((sqrt(h[1]/n0) * (matrix(K_h(x0,X[i,1],h[1]),ncol=1) /pjx +  Qest %*%Khh))^2)}
  
  
  registerDoMC(cores=3)
  tmp2=foreach(i = 1:n0, .combine='+')  %dopar% serror_fun(i)
  
  obj1=est_f_debiased-qnorm(0.975)*sqrt(tmp2/(n0*h[1]))  
  obj2=est_f_debiased+qnorm(0.975)*sqrt(tmp2/(n0*h[1]))  
  
  results_set=cbind(results_set, sqrt(n0*h[1]) * (est_f_debiased-f_eval[,j_true])/sqrt(tmp2))
}









