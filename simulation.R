source("generateModel.R")
source("ricf_int2.R")
source("ricf_int.R")
source("ricf_dg.R")

#library(igraph)


# v: graph size
# n: multiple of sample size to graph size
# l: length of the unique cycle
# d: probability of filling an edge between two nodes
replicates <- 1000
V <- c(5, 10, 20)
N <- c(2.5, 5, 10, 20)
K <- c(0, 3, 4, 8)
D <- c(0.2)

v <- V[1]
n <- N[3]
l <- K[2]
d <- D[1]

rss <- function(r, models, List){
  rss <- sum((t(models[[r]]$B.true)-List[[r]]$Lambdahat)^2) + sum((diag(models[[r]]$Omega.true)-List[[r]]$Omegahat)^2)
  num <- sum(models[[r]]$B) + nrow(models[[r]]$Omega)
  return (rss / num)
}

d_llh <- function(r, models, List){
  #sigma <- solve(diag(v) - models[[r]]$B.true) %*% models[[r]]$Omega.true %*% t(solve((diag(v)-models[[r]]$B.true)))
  #K <- t(diag(v) - models[[r]]$B.true) %*% solve(models[[r]]$Omega.true) %*% (diag(v) - models[[r]]$B.true)
  K = (diag(v)-List[[r]]$Lambdahat) %*% (1 / List[[r]]$Omegahat * t(diag(v)-List[[r]]$Lambdahat))
  tr <- sum( c(K) * c(t(models[[r]]$Sigmas[[1]])) )
  llh <- - sum(log(List[[r]]$Omegahat)) + log(det(diag(v)-List[[r]]$Lambdahat) ^ 2) - tr
  llh0 <- - sum(log(diag(models[[r]]$Omega.true))) + log(det(diag(v)-models[[r]]$B.true) ^ 2) - v
  return (llh0 - llh)
}

run_vnl <- function(v, n, l, d=0.2, seed = 1, repl = 1000){
  set.seed(seed)
  graphs <- list()
  models <- list()
  
  for (r in 1:repl){
    graphs[[r]] <- generateGraph(v=v, l=l, d=d, b=0)
  }
  
  for (r in 1:repl){
    targets <- list(numeric(0))
    target.length <- c(round(n*v))
    for (k in 1:sample(3,1)){
      targets[[k+1]] <- sample(v, sample(2:3,1))
      target.length[k+1] <- max(sample(round(n*v/2), 1) + round(n*v/2),v+1)
    }
    if (r %% 100 == 0)
      print(r)
    #models[[r]] <- generateModel(v=v, l=l, d=d, b=0, targets = targets, target.length = target.length)
    models[[r]] <- generateData(B=graphs[[r]]$B, Omega=graphs[[r]]$Omega, targets = targets, target.length = target.length)
  }
  
  res_int <- list()
  res_dg_agg <- list()
  
  t0 <- proc.time()[3]
  for (r in 1:repl){
    # preprocess, to simplify variable names
    # targets, target.length, data, empty matrix/vector for recording estimates
    # count_L, count_O: sum of valid data rows for each position in zero pattern matrix/vector
    targets <- models[[r]]$targets
    target.length <- models[[r]]$target.length
    data <- t(models[[r]]$Y)
    Lhat <- count_L <- matrix(0,v,v)
    Omegahat <- count_O <- rep(0,v)
    
    # the zero pattern of L,O, observation case
    # allconverged: whether all ricf in the aggregation converged
    L0 <- t(models[[r]]$B)
    Omega0 <- rep(1,v)
    allconverged <- T
    for (k in 1:length(targets)){
      # zero pattern of intervened graphs
      L <- L0
      L[, targets[[k]]] <- 0
      Omega <- Omega0
      Omega[targets[[k]]] <- 0
      
      # find the indices of rows in data for intervention target k
      # assume that length(ind) > 0 always holds
      posind <- cumsum(c(0,target.length))
      ind <- (posind[k]+1): posind[k+1]
      m <- length(ind)
      S <- cov(data[ind, ])
      #print(n)
      
      # once ricf performing, weighted sum of Lambdahat, Omegahat * zero pattern of intervened graphs
      # update allconverged
      tmp <- ricf_R_dg(L=L, S=S)
      Lhat <- Lhat + m * tmp$Lambdahat * L
      Omegahat <- Omegahat + m * tmp$Omegahat * Omega
      count_L <- count_L + m * L
      count_O <- count_O + m * Omega
      allconverged <- tmp$converged & allconverged
    }
    
    # weighted average of K ricf results
    # set the zeros in counts to be 1
    count_L <- pmax(count_L, 1)
    count_O <- pmax(count_O, 1)
    res_dg_agg[[r]] <- list()
    res_dg_agg[[r]]$Lambdahat <- Lhat / count_L
    res_dg_agg[[r]]$Omegahat <- Omegahat / count_O
    # vector * matrix is a row transformation
    res_dg_agg[[r]]$Sigmahat <-  t(solve((diag(v)-res_dg_agg[[r]]$Lambdahat))) %*% (res_dg_agg[[r]]$Omegahat * solve((diag(v)-res_dg_agg[[r]]$Lambdahat)))
    res_dg_agg[[r]]$allconverged <- allconverged
  }
  t1 <- (proc.time()[3] - t0) #/ replicates
  
  t0 <- proc.time()[3]
  for (r in 1:repl){
    res_int[[r]] <- ricf_int_(L=t(models[[r]]$B), data=t(models[[r]]$Y), targets=models[[r]]$targets, target.length=models[[r]]$target.length)
  }
  t2 <- (proc.time()[3] - t0) #/ replicates
  
  conv_mle <- sum(sapply(1:repl, function (r) res_int[[r]]$converged == T))
  conv_agg <- sum(sapply(1:repl, function (r) res_dg_agg[[r]]$allconverged == T))
  rss_mle <- mean(sapply(1:repl, rss, models=models, List=res_int))
  rss_agg <- mean(sapply(1:repl, rss, models=models, List=res_dg_agg))
  llh_mle <- mean(sapply(1:repl, d_llh, models=models, List=res_int))
  llh_agg <- mean(sapply(1:repl, d_llh, models=models, List=res_dg_agg))
  
  res <- matrix(c(conv_agg, conv_mle, rss_agg, rss_mle,
                  llh_agg, llh_mle, t1, t2), 2, 4)
  return (res)
}

for (v in V){
  for (n in N){
    for (l in K){
      run_vnl(v, n, l)
    }
  }
}



tmp=(sapply(1:replicates, rss, List=res_int2))
ind = which(tmp>10)
for (r in ind){
  print(models[[r]]$targets)
  cat("-----------------")
  plotGraph(models[[r]]$B)
  Sys.sleep(10)
}
