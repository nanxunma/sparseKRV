library(mvtnorm)
library(Rcpp)
cppFunction('
  double ffunction(NumericMatrix Z, NumericVector c, double rho) {
    int n = Z.nrow();
    double f = 0;
    double inner = 0;
    for (int i = 0; i < n; i++){
      for (int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        f = f + exp(-2.0/rho*inner);
      }
    }
    return f;
  }

  ')


cppFunction('
  double hfunction(NumericMatrix Z, NumericVector c, NumericMatrix Ky, double rho){
    int n = Z.nrow();
    double h = 0.0;
    double inner = 0.0;
    for (int i = 0; i < n; i++){
      for (int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        h = h + Ky(i, ip) * exp(-1.0 / rho * inner);
      }
    }

    return h;

  }
  ')


cppFunction('
  double fgrad(int k, NumericMatrix Z, NumericVector c, double rho){
    int n = Z.nrow();
    double f = 0.0;
    double inner = 0.0;
    for (int i = 0; i < n; i++){
      for (int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        f = f + exp(-2 / rho * inner ) * (-2.0 / rho) * pow(Z(i,k)-Z(ip,k), 2.0);
      }
    }
    return f;

  }
  ')


cppFunction('
  NumericVector fgradall(NumericMatrix Z, NumericVector c, double rho){
    int n = Z.nrow();
    int p = Z.ncol();
    NumericVector f (p);
    double inner = 0;
    for(int i = 0; i < n; i++){
      for(int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        f = f + exp(-2/rho*inner) * (-2/rho) * pow(Z(i,_)-Z(ip,_),2);
      }
    }
    return f;
  }')



cppFunction('
  double hgrad(int k, NumericMatrix Z, NumericVector c, NumericMatrix Ky, double rho){
    int n = Z.nrow();
    double h = 0;
    double inner = 0;
    for(int i = 0; i < n; i++){
      for(int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        h = h + Ky(i, ip) * exp(-1.0/rho * inner ) * (-1/rho) * pow(Z(i,k)-Z(ip,k), 2);
      }
    }
    return h;
  }')


cppFunction('
  NumericVector hgradall(NumericMatrix Z, NumericVector c, NumericMatrix Ky, double rho){
    int n = Z.nrow();
    int p = Z.ncol();
    NumericVector h (p);
    double inner = 0;
    for(int i = 0; i < n; i++){
      for(int ip = 0; ip < n; ip++){
        inner = sum(c * pow(Z(i,_)-Z(ip,_), 2.0));
        h = h + Ky(i, ip) * exp(-1.0/rho * inner ) * (-1/rho) * pow(Z(i,_)-Z(ip,_), 2);
      }
    }
    return h;
  }')



Lgrad <- function(x, k, c, Z, Ky, rho=1, w, lambda, mu, z, M){
  c[k] = x
  L = fgrad(k, Z, c, rho) + mu*hgrad(k, Z, c, Ky, rho) + w * (hfunction(Z, c, Ky, rho) - M) * hgrad(k, Z, c, Ky, rho) + lambda[k] + w*c[k] - w*z[k]
  return(L)
}


Lfunction <- function(x, k, c, Z, Ky, rho=1, w, lambda, mu, z, M){
  c[k] = x
  f = ffunction(Z,c,rho)
  h = hfunction(Z,c,Ky,rho)
  L = f + mu*h + w/2*(h-M)^2 + t(lambda)%*%c + w/2 * t(c-z)%*%(c-z)
  return(L)
}

soft <- function(a,k){
      rst <- sign(a)*pmax(0, abs(a)-k)
      return(rst)
}


admm.gaussian <- function(Y, Z, Ky, rho, s, w=1, w.increment = 1, maxiter = 200, tol=1e-3, cinit=NULL, zinit=NULL){
  n = dim(Y)[1]
  p = dim(Z)[2]
  r = dim(Y)[2]
  iter = 0
  diff_value = diff2 = 10
  ## initial value
  #c = abs(rnorm(p,mean=1, sd = 1))
  if(is.null(cinit)) {c = rnorm(p,0,0.01)} else{
    c = cinit
  }
    if(is.null(zinit)) {z = rnorm(p)} else{
    z = zinit
  }
  mu = rnorm(1)
  lambda = rnorm(p)
  I = diag(rep(1,p))
  #Ky = tcrossprod(Y)
  M = n*sqrt(sum(diag(crossprod(Ky))))
  while((iter==0) || (iter<maxiter) && (diff_value > tol | diff2 > tol))
  {
    cat(iter)
    # update c
    c.prev = c
    for(k in 1:p){
      #c[k] = optim(0.1, function(x) Lfunction(x,k,c,Z,Ky,rho,w,lambda,mu,z,M), method = "Brent",lower=0, upper=M)$par[1]
      #cat(k)
      #c[k] = optim(c.prev[k], function(x) Lfunction(x,k,c,Z,Ky,rho,w,lambda,mu,z,M), method = "Brent",lower=0, upper=M)$par[1]
      c[k] = optim(c.prev[k], function(x) Lfunction(x,k,c,Z,Ky,rho,w,lambda,mu,z,M), method = "L-BFGS-B",lower=0, upper=M)$par[1]
      #c[k] = optim(c.prev[k], function(x) Lfunction(x,k,c,Z,Ky,rho,w,lambda,mu,z,M), method = "BFGS",lower=0, upper=M)$par[1]
      #c[k] = optim(0.1, function(x) abs(Lgrad(x,k,c,Z,Ky,rho,w,lambda,mu,z, M)), method="Brent",lower=0, upper=M)$par[1]
      #c[k] = optim(c.prev[k], function(x) abs(Lgrad(x,k,c,Z,Ky,rho,w,lambda,mu,z, M)), method="BFGS",lower=0, upper=M)$par[1]
      if(is.na(c[k])){c[k]<-c.prev[k]}
    }
    # update z
    z.prev = z
    z = soft(c+lambda/w, s/w)
    # update mu
    mu.prev = mu
    hh = hfunction(Z, c, Ky, rho)
    mu = mu + w*(hh - M)
    #cat("hfunction is", hh)
    # update lambda
    lambda.prev = lambda
    lambda = lambda + w*(c-z)
    # update iteration and difference
    iter = iter + 1
    #diff_value = sum((c-z)^2) + abs(t(c)%*%g-1)
    diff_value = sum(abs(c-c.prev)) / (sum(abs(c.prev))+0.0001)
    diff2 = sum(abs(c-z))
    w = w*w.increment
  }
  z[abs(z)<1e-10] <- 0
  diff1 = hfunction(Z, c, Ky, rho) - M
  out <- list(c=c, z=z, diff1 = diff1, diff2=diff2, diff0=diff_value, iters=iter, sparse=sum(z!=0))
  return(out)
}

beta.matrix <- function(r, p, d, mm, mean = 1){
  B <- matrix(0, nrow = r, ncol = p)
  for(i in 1:d){
    index = 1:mm
    #B[index,i] <- rnorm(mm, mean = mean, sd=0.2)
    B[index,i] <- mean
  }
  if(d < r){
    for(j in (d+1):p){
      #B[,j] <- rnorm(r, mean = 0.01, sd = 0.01)
      B[,j] <- rep(0,r)
    }}
  return(B)
}


beta.matrix2 <- function(r, p, d, mm, mean = 1){
  B <- matrix(0, nrow = r, ncol = p)
  for(i in 1:d){
    index = sample(1:r, mm, replace = F)
    B[index,i] <- 1
  }
  if(d < r){
    for(j in (d+1):p){
      #B[,j] <- rnorm(r, mean = 0.01, sd = 0.01)
      B[,j] <- rep(0,r)
    }}
  return(B)
}

KRV.function <- function(Ky, Z, c) {
  Z0 = Z[,c!=0]
  Kz = tcrossprod(Z0)
  #krv = sum(diag(Ky %*% Kz)) / (sqrt(sum(diag(Ky %*% Ky)))*sqrt(sum(diag(Kz %*% Kz))))
  krv = sum(diag(crossprod(Ky,Kz))) / sqrt(sum(diag(crossprod(Ky)))*sum(diag(crossprod(Kz))))
  return(krv)
}

KRV.linear <- function(Ky,Z,c){
  n = dim(Z)[1]
  Kz = matrix(NA,n,n)
  for(i in 1:n){
    for(ip in 1:n){
      Kz[i,ip] = sum(c*Z[i,]*Z[ip,])
    }
  }
  krv = sum(diag(crossprod(Ky,Kz))) / sqrt(sum(diag(crossprod(Ky)))*sum(diag(crossprod(Kz))))
  return(krv)
}


rm.zero <- function(K){
  K[is.na(K)] <- 0
  return(K)
}


Kz.gaussian <- function(Z, c, rho){
  n = dim(Z)[1]
  Kz = matrix(NA,n,n)
  for(i in 1:n){
    for(ip in 1:n){
      Kz[i,ip] = exp(sum(-c*(Z[i,]-Z[ip,])^2/rho))
    }
  }
  return(Kz)
}





#library(MiSPU)
library(GUniFrac)
library(dirmult)
library(PMA)
library(cluster)
#library(phangorn)

  data(throat.otu.tab)
  
  otutab= throat.otu.tab
  
  ## DirMultOutput data
  dd = dirmult(otutab)
  data(throat.tree)





do.one <- function(n, p, r, d, s0, mm, mean=5, var=1){
  

  S = simPop(J=n,  n= 1000, pi= dd$pi, theta= dd$theta)
  Y = S$data
  colnames(Y) <- names(S$pi)
  err <- matrix(rnorm(n*p, 0, 1), ncol = p, nrow = n)
  
#### bray
  Ky1 = as.matrix(vegdist(Y,method="bray"))
  Ky1 = rm.zero(Ky1)
  #Ky1 = psd.kernel(Ky1)
  #r = dim(Y)[2]
  OTU.mean = apply(Y,2,mean)
  YJ = Y[,order(OTU.mean,decreasing = T)]
  std = apply(YJ,2,mean)
  Y2 = YJ %*% diag(1/std)
  
  #Y = rmvnorm(n=n, mean = c(rep(mean,s0),rep(0,r-s0)), sigma = diag(rep(1,r)))
  beta = beta.matrix(r,p,d,mm=mm,mean = 1)
 
  Z <- YJ %*% beta
  Z[,1:d] <- scale(Z[,1:d])
  Z <- Z + err
  #  Z <- rmvnorm(n=n, mean = c(apply(Y,2,mean),rep(0,p-r)),sigma = diag(rep(var,p)))
  In = diag(rep(1,n))
  Ky = (In - 1/n) %*% Ky1 %*% Ky1 %*% (In-1/n)

  M = n*sqrt(sum(diag(crossprod(Ky))))
  
 s.seq = 10^seq(-2,6,by=0.02)
 rst.krv = rep(NA,length(s.seq))
 rst.krv = matrix(NA, ncol=5, nrow=length(s.seq))
 folds <- cut(seq(1,n),breaks=5,labels=FALSE)
 for(i in 1:5){
   testIndexes <- which(folds==i,arr.ind=TRUE)
   testY <- Y[testIndexes, ]
   trainY <- Y[-testIndexes, ]
   testZ <- Z[testIndexes, ]
   trainZ <- Z[-testIndexes, ]
   M.5fold[i] = n*sqrt(sum(diag(crossprod(tcrossprod(trainY)))))
 }
 
 krv.all = rep(NA, length(s.seq))
 for(kk in 1:length(s.seq)){
   tun = s.seq[kk]
   for(i in 1:5){
     testIndexes <- which(folds==i,arr.ind=TRUE)
     testY <- Y[testIndexes, ]
     trainY <- Y[-testIndexes, ]
     testZ <- Z[testIndexes, ]
     trainZ <- Z[-testIndexes, ]
     n1 = dim(testZ)[1]
     n2 = dim(trainZ)[1]
          sim.Ky = as.matrix(vegdist(testY,method="bray"))
     Ky.test = rm.zero(sim.Ky)
     Ky.test = (In - 1/n1) %*% Ky.test %*% Ky.test %*% (In-1/n1)

          sim.Ky =as.matrix(vegdist(trainY,method="bray"))
     Ky.train = rm.zero(sim.Ky)
     Ky.train = (In - 1/n2) %*% Ky.train %*% Ky.train %*% (In-1/n2)

     #fit.train = admm.linear(trainY,trainZ,s=n*tun, w=1,w.increment=1.1,maxiter = 1000, tol = 1e-4)
     fit.train = admm.gaussian(Y=trainY,Z=trainZ,Ky=Ky.train,s=n2*tun, M=M.5fold[i],w=1,w.increment = 1.1,maxiter = 1000,tol = 1e-4)
     fit.z = fit.train$z
     rst.krv[kk,i] <- KRV.function(Ky.test, testZ, fit.z)
   }

 }
 krv.mean <- apply(rst.krv,1,mean)
 krv.local <- apply(rst.krv,2, max)

 s.optim = max(s.seq[which(krv.mean==max(krv.mean))]) + 2*sd(krv.local)
 #s.optim = max(chosen) + 2*sd(krv.local)
 
  print(s.optim)
  out.bray = admm.gaussian(Y=Y,Z=Z,Ky=Ky,rho=rho,s=n*s.optim, w=1,w.increment=1.1,maxiter = 1000, tol = 1e-4)

### bray end


  YC = Y[,OTU.mean !=0]
    perm.out = CCA.permute(YC,Z, typex="standard",typez="standard", nperms=100)
    out <- CCA(YC,Z,typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init)

   s.seq = 10^seq(-3,0,by=0.02)
   rst.cca = matrix(NA, ncol=5, nrow=length(s.seq))
   folds <- cut(seq(1,n),breaks=5,labels=FALSE)
   for(kk in 1:length(s.seq)){
     tun = s.seq[kk]
     for(i in 1:5){
       testIndexes <- which(folds==i,arr.ind=TRUE)
       testY <- YC[testIndexes, ]
       trainY <- YC[-testIndexes, ]
       testZ <- Z[testIndexes, ]
       trainZ <- Z[-testIndexes, ]
       tt.index = ((apply(testY,2,sum)!=0) +( apply(trainY,2,sum)!=0)) == 2
       testY = testY[,tt.index]
       trainY = trainY[,tt.index]
      cca.cv = CCA(trainY,trainZ, typex="standard",typez="standard",penaltyx=NULL,penaltyz=tun)
     u1 = cca.cv$u
     v1 = cca.cv$v  
 cca.stat = t(u1) %*% t(testY) %*% testZ %*% v1
 rst.cca[kk,i] <- cca.stat[1,1]
}
}
final.cca <- apply(rst.cca,1,mean)
tun.cca.optim <- max(s.seq[which(final.cca==min(final.cca))]) 
cca.cv.f = CCA(trainY,trainZ, typex="standard",typez="standard",penaltyx=NULL,penaltyz=tun.cca.optim)

    return(c(as.vector(out.bray$z), as.vector(cca.cv.f$v)))
}


result = replicate(nSim, do.one(n=n, p=p, r=r, d=d, s0 = s0, mm=mm, mean = mean))
name = paste("bray-gaussian.most10.n",n,".p",p,".r",r,".d",d,".rho",rho,".csv", sep = "")
write.csv(result, name)
