library(mvtnorm)
soft <- function(a,k){
  rst <- sign(a)*pmax(0, abs(a)-k)
  return(rst)
}
# generate a psd kernel matrix
Posdef <- function (n, ev = runif(n, 0, 10)) 
{ 
  Z <- matrix(ncol=n, rnorm(n^2)) 
  decomp <- qr(Z) 
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp) 
  d <- diag(R) 
  ph <- d / abs(d) 
  O <- Q %*% diag(ph) 
  Z <- t(O) %*% diag(ev) %*% O 
  return(Z) 
} 


B.matrix <- function(Z){
  n = dim(Z)[1]
  p = dim(Z)[2]
  B = matrix(NA,p,p)
  for(j in 1:p){
    for(k in 1:p){
      M = matrix(0,n,n)
      for(i in 1:n){
        for(ip in 1:n){
          M[i,ip] = Z[i,j]*Z[ip,j]*Z[i,k]*Z[ip,k]
        }
      }
      B[j,k] = sum(M)
    }
  }
  return(B)
}


g.vector <- function(Z, Ky){
  n = dim(Z)[1]
  p = dim(Z)[2]
  g = rep(NA, p)
  for(j in 1:p){
    g[j] = 0
    for(i in 1:n){
      for(ip in 1:n){
        g[j] = g[j] + Ky[i,ip]*Z[i,j]*Z[ip,j]
      }
    }
  }
  return(g)
}

admm.linear.inside = function(B, g, s, M, w, w.increment, maxiter, tol)
{
  p = length(g)
  iter = 0
  diff_value = diff2 = 10
  ## initial value
  c = rnorm(p)
  z = rnorm(p)
  lambda = rnorm(p)
  mu = rnorm(1)
  I = diag(rep(1,p))
  while((iter==0) || (iter<maxiter) && (diff_value > tol | diff2 > tol))
  {
    # update c
    c.prev = c
    c = solve(2*B + w*(g%*%t(g) + I)) %*% (as.numeric(w*M - mu)*g + w*z - lambda)
    # update z
    z.prev = z
    z = soft(c+lambda/w, s/w)
    # update mu
    mu.prev = mu
    mu = mu + w*(t(c)%*%g - M)
    # update lambda
    lambda.prev = lambda
    lambda = lambda + w*(c-z)
    # update iteration and difference
    iter = iter + 1
    #diff_value = sum((c-z)^2) + abs(t(c)%*%g-1)
    diff_value = sum(abs(c-c.prev)) / sum(abs(c.prev))
    diff2 = sum(abs(c-z))
    w = w*w.increment
  }
  diff1 = sum(abs(c-z))
  diff2 = t(c)%*%g -M
  out <- list(c=c, z=z, diff1 = diff1, diff2=diff2, diff0=diff_value, iters=iter, sparse=sum(z!=0))
  return(out)
}




admm.linear <- function(Y, Z, Ky, s, w=1, w.increment = 1, maxiter = 1000, tol=1e-4){
  n = dim(Y)[1]
  p = dim(Z)[2]
  r = dim(Z)[2]
  #Ky = tcrossprod(Y)
  M = sqrt(sum(diag(crossprod(Ky))))
  B = B.matrix(Z)
  g = g.vector(Z, Ky)
  out <- admm.linear.inside(B,g,s=s, w=w, M=M,w.increment=w.increment,maxiter = maxiter, tol = tol)
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

psd.kernel <- function(K){
  eig <- eigen(K)
  K2 = eig$vectors %*% diag(abs(eig$values)) %*% t(eig$vectors)
  return(K2)
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
  

  sim.Ky = GUniFrac(Y, throat.tree, alpha=c(0, 0.5, 1))$unifracs
  Ky1 = as.matrix(vegdist(Y,method="bray"))
  Ky1 = rm.zero(Ky1)
  #Ky1 = psd.kernel(Ky1)
  #r = dim(Y)[2]
  OTU.mean = apply(Y,2,mean)
  YJ = Y[,order(OTU.mean,decreasing = T)]
  std = apply(YJ,2,mean)
  Y2 = YJ %*% diag(1/std)

##clustering
  cluster = pam(t(Y),20, metric="euclidean")
  summary(cluster)
  cluster$clusinfo
  id = cluster$clustering

  member = matrix(0, ncol=2, nrow=20)
  member[,1] = 1:20
  for(i in 1:20){
    class = Y[,id == i]
    cat("class",i, "is", sum(class!=0)/sum(Y!=0),"\n")
    member[i,2] = sum(class!=0)/sum(Y!=0)
  }
  order.m = order(member[,2],decreasing=T)
  kclass = order.m[2]
  Y.toadd = Y[,id == kclass]
  d.toadd = dim(Y.toadd)[2]

  YJ[,1:d.toadd] = Y.toadd
  #Y = rmvnorm(n=n, mean = c(rep(mean,s0),rep(0,r-s0)), sigma = diag(rep(1,r)))
  beta = beta.matrix(r,p,d,mm=d.toadd,mean = 1)
  err <- matrix(rnorm(n*p, 0, 1), ncol = p, nrow = n)
  Z <- YJ %*% beta
  Z[,1:d] <- scale(Z[,1:d])
  Z <- Z + err
  #  Z <- rmvnorm(n=n, mean = c(apply(Y,2,mean),rep(0,p-r)),sigma = diag(rep(var,p)))
  In = diag(rep(1,n))
  Ky = (In - 1/n) %*% Ky1 %*% (In-1/n)
  B = B.matrix(Z)
  g = g.vector(Z, Ky)
  M = n*sqrt(sum(diag(crossprod(Ky))))

 s.seq = 10^seq(-2,6,by=0.02)
 rst.krv = rep(NA,length(s.seq))
 rst.krv = matrix(NA, ncol=5, nrow=length(s.seq))
 B.5fold <- array(NA, dim=c(5,p,p))
 g.5fold <- array(NA, dim=c(5,p))
 M.5fold <- rep(NA,5)
 folds <- cut(seq(1,n),breaks=5,labels=FALSE)
 for(i in 1:5){
   testIndexes <- which(folds==i,arr.ind=TRUE)
   testY <- Y[testIndexes, ]
   trainY <- Y[-testIndexes, ]
   testZ <- Z[testIndexes, ]
   trainZ <- Z[-testIndexes, ]
   B.5fold[i,,] <- B.matrix(trainZ)
   sim.Ky = GUniFrac(trainY, throat.tree, alpha=c(0, 0.5, 1))$unifracs
   Ky1 = sim.Ky[,,"d_UW"]
   g.5fold[i,] <- g.vector(trainZ, Ky1)
   
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
     

     fit.train = admm.linear.inside(B=B.5fold[i,,],g=g.5fold[i,],s=n/4*5*tun, M=M.5fold[i],w=1,w.increment = 1.1,maxiter = 1000,tol = 1e-4)
     fit.z = fit.train$z
     sim.Ky = GUniFrac(testY, throat.tree, alpha=c(0, 0.5, 1))$unifracs
     Ky.test= sim.Ky[,,"d_UW"]
     rst.krv[kk,i] <- KRV.linear(Ky.test, testZ, fit.z)
 
 
   }

   fit.all = admm.linear.inside(B=B,g=g,s=n*tun,M=M,w=1,w.increment = 1.1,maxiter = 1000,tol=1e-4)

 }
 krv.mean <- apply(rst.krv,1,mean)
 krv.local <- apply(rst.krv,2, max)

 
 s.optim = max(s.seq[which(krv.mean==max(krv.mean))]) + 2*sd(krv.local)
 #s.optim = max(chosen) + 2*sd(krv.local)
 print(s.optim)
 out2 = admm.linear.inside(B=B,g=g,s=n*s.optim,M=M, w=1,w.increment=1.1,maxiter = 1000, tol = 1e-4)
 # # #### block 2 end
 #t2 = CCA(YJ,Z, typex="standard",typez="standard")
  YC = Y[,OTU.mean !=0]
    perm.out = CCA.permute(YC,Z, typex="standard",typez="standard", nperms=100)
    out <- CCA(YC,Z,typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init)
    return(c(as.vector(out2$z), as.vector(out$v)))
}


result = replicate(nSim, do.one(n=n, p=p, r=r, d=d, s0 = s0, mm=mm, mean = mean))
name = paste("linear.bray.n",n,".p",p,".r",r,".d",d,".s0",s0,".mm",mm,".mean",mean,".csv", sep = "")
write.csv(result, name)
