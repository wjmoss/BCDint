source('generateModel.R'); source('ricf_int.R'); source('ricf_dg.R')
set.seed(19)
# Observation-only one-variable closed form.
x <- matrix(c(1,3,5,7),4,1)
f <- ricf_int_(L=matrix(0,1,1),data=x,maxiter=1)
stopifnot(f$converged, f$Bhat[1,1]==1,
          identical(as.numeric(f$Omegahat),var(as.numeric(x))))
# Intervened rows must not contribute to shared observational noise variance.
x2 <- rbind(x,matrix(c(100,300,500),3,1))
g <- ricf_int_(L=matrix(0,1,1),data=x2,
              targets=list(numeric(0),1),target.length=c(4,3))
stopifnot(identical(as.numeric(g$Omegahat),var(as.numeric(x))))
# Generation with explicit n must not read a global sample count.
n <- 999;v <- 999
B <- matrix(0,3,3);B[2,1]<-B[3,2]<-1
m <- generateData(B,diag(3),n=20)
stopifnot(identical(dim(m$Y),c(3L,20L)),identical(m$targets,list(numeric(0))))
m2 <- generateData(B,diag(3),targets=list(numeric(0),2),n=20)
stopifnot(identical(dim(m2$Y),c(3L,40L)),all(m2$target.length==20))
tryCatch({generateData(B,diag(3));stop('Expected count validation error')},
 error=function(e) stopifnot(grepl('Supply target.length',conditionMessage(e))))
# Independent covariance-inverse likelihood evaluation.
f2 <- ricf_int(L=t(B),data=t(m2$Y),targets=m2$targets,
               target.length=m2$target.length,restarts=2,tol=1e-9)
reference<-0;off<-c(0,cumsum(m2$target.length))
for(k in seq_along(m2$targets)) {
 y<-t(m2$Y)[seq.int(off[k]+1,off[k+1]),,drop=FALSE]
 L<-f2$Lambdahat;L[,m2$targets[[k]]]<-0
 O<-f2$Omegahat;O[m2$targets[[k]]]<-diag(cov(y))[m2$targets[[k]]]
 inv<-solve(diag(3)-L);Sigma<-t(inv)%*%diag(O)%*%inv
 reference<-reference+nrow(y)*(-as.numeric(determinant(Sigma)$modulus)-sum(diag(solve(Sigma,cov(y)))))
}
stopifnot(abs(reference-llh_int(f2,t(m2$Y),m2$targets,m2$target.length))<1e-8,
 f2$joint_score==max(f2$restart_scores))
h<-ricf_int(B=B,data=t(m$Y),restarts=0,maxiter=1)
stopifnot(h$iterations==1)
cat('Passed single-variable, n/count validation, joint score, restart and control checks.\n')
# The exact zero-mean likelihood must agree with a direct covariance evaluation.
u <- ricf_int(L=t(B),data=t(m2$Y),targets=m2$targets,
              target.length=m2$target.length,restarts=2,covariance="ml")
reference_ml <- 0
for(k in seq_along(m2$targets)) {
 y<-t(m2$Y)[seq.int(off[k]+1,off[k+1]),,drop=FALSE]
 L<-u$Lambdahat;L[,m2$targets[[k]]]<-0
 S<-crossprod(y)/nrow(y)
 O<-u$Omegahat;O[m2$targets[[k]]]<-diag(S)[m2$targets[[k]]]
 inv<-solve(diag(3)-L);Sigma<-t(inv)%*%diag(O)%*%inv
 reference_ml<-reference_ml+nrow(y)*(-as.numeric(determinant(Sigma)$modulus)-sum(diag(solve(Sigma,S))))
}
stopifnot(abs(reference_ml-u$joint_score)<1e-8,u$covariance=="ml")
cat('Passed exact Gaussian likelihood check.\n')
