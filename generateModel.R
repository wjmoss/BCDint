#' Generate Data from Structural Equation Model
#'
#'
#' Samples a random Graph and Data
#'
#'
#' @param v number of nodes
#' @param n number of samples
#' @param l the size of the largest cycle
#' @param d the
#' @return \item{sigmaHat}{estimated covariance matrix at convergence}
#'    \item{bHat}{estimated B matrix (edge weights for directed edges) at converegence}
#'    \item{omegaHat}{estimated Omega (edge weights for bi-directed edges) at convergence}
#'    \item{iterations}{number of iterations until convergence. a single iteration is considered
#'    a pass through all nodes}
#'    \item{converged}{boolean on whether or not the algorithm converged before the max iterations}
generateModel <- function(v, n=NULL, l, d, b=0, scc_constraint=T, acy_constraint=F, det_constraint=0.01,
                          targets=NULL, target.length=NULL, errorDist = "gauss", errorInt="gauss"){
  # Setup B and Omega matrices
  B <- matrix(0, nrow = v, ncol = v)
  Omega <- diag(rep(1, v))

  # Generate cycle from 1 -> 2 -> ... K -> 1
    if(l > 1){
      for(i in 1:(l-1))
      {
        B[i+1, i] <- 1
      }
      B[1, l] <- 1
    }

    # fill in remaining edges
    for(j in 1:(v-1))
    {
      for(i in max(j+1, l+1):v)
      {
        if((B[i,j] != 1) & (B[j,i]!= 1)){
          U <- runif(1)
          if(U < d){
            B[i, j] <- 1
          } else if(U < (b + d))
          {
            Omega[i,j] <- Omega[j,i] <- 1
          }
        }
      }
    }


  # reorder the vertices
  reorder <- sample(v)
  B <- B[reorder, reorder]
  Omega <- Omega[reorder, reorder]

  repeat{
    # Sample edge weights as given in the paper
    B.true <- matrix(runif(v^2,0.3,1), nrow = v) * B
    B.true <- matrix(rbinom(v^2,1,0.5)*2-1, nrow = v) * B.true
    
    if (det(diag(v) - B.true) > det_constraint)
      break
  }
  
  Omega.true <- matrix(runif(v^2,-1,1), nrow = v)
  Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)]
  Omega.true <- Omega.true * Omega

  # ensure Omega.true is PD by making it diagonally dominant
  #for(i in 1:v){
    #Omega.true[i, i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
  #}
  diag(Omega.true) <- runif(v, 0.3, 1)

  data <- list()
  Y <- NULL
  Sigmas <- list()
  if (is.null(targets) || is.null(target.length)){
    targets = list(numeric(0))
    target.length = (n)
  }

  for (k in 1:length(target.length)){
    tar <- targets[[k]]
    n <- target.length[k]
    B.k <- B.true
    B.k[tar, ] <- 0
    Omega.k <- Omega.true
    Omega.k[tar, ] <- Omega.k[, tar] <- 0
    for (t in tar){
      for (i in 1:v){
         if (Omega[t, i] != 0 && !any(tar == i)){
          Omega.k[i, i] <- Omega.k[i, i] + sqrt(abs(Omega.true[t, i]*sqrt(Omega.true[i, i]*Omega.true[t, t])))
        }
      }
    }
    diag(Omega.k)[tar] <- diag(Omega.true)[tar]

    # Sample data from multivariate normal and make mean 0
    if(errorDist == "gamma")
    {
      errors <- matrix(rgamma(v * n, shape = 1, scale = 1), nrow = v, ncol = n)
    } else if(errorDist == "poisson")
    {
      pois.mu <- 3
      errors <- matrix(rpois(v * n, lambda = pois.mu), nrow = v, ncol = n)
      errors <- (errors - pois.mu)/sqrt(pois.mu)
    } else if(errorDist == "t"){
      t.df <- 5
      errors <- matrix(rt(v * n, df = t.df) * (t.df - 2)/t.df, nrow = v, ncol = n)
    } else if(errorDist == "lognormal"){
      ln.var <- 1
      ln.mu <- 0
      errors <- matrix(rlnorm(v * n, meanlog = ln.mu, sdlog = sqrt(ln.var)), nrow = v, ncol = n)
      errors <- (errors - exp(ln.mu + ln.var/2)) / sqrt((exp(ln.var) - 1) * exp(2*ln.mu + ln.var))
    } else if(errorDist == "gauss") {
      errors <- matrix(rnorm(v * n, mean = 0, sd = 1), nrow = v, ncol = n)
      #Y <- t(MASS::mvrnorm(n = n, mu = rep(0, v), Sigma = sigma))
    }

    #cov(t(errors))
    temp <- solve(diag(rep(1,v)) - B.k)
    sigma.k <- temp %*% Omega.k %*% t(temp)
    errors <- t(chol(Omega.k)) %*% errors

    # interventional errors distribution
    if(errorInt == "gamma")
    {
      errors <- matrix(rgamma(length(tar) * n, shape = 1, scale = 1), nrow = length(tar), ncol = n)
    } else if(errorInt == "poisson")
    {
      pois.mu <- 3
      errors[tar, ] <- matrix(rpois(length(tar) * n, lambda = pois.mu), nrow = length(tar), ncol = n)
      errors <- (errors - pois.mu)/sqrt(pois.mu)
    } else if(errorInt == "t"){
      t.df <- 5
      errors[tar, ] <- matrix(rt(length(tar) * n, df = t.df) * (t.df - 2)/t.df, nrow = length(tar), ncol = n)
    } else if(errorInt == "lognormal"){
      ln.var <- 1
      ln.mu <- 0
      errors[tar, ] <- matrix(rlnorm(length(tar) * n, meanlog = ln.mu, sdlog = sqrt(ln.var)), nrow = length(tar), ncol = n)
      errors[tar, ] <- (errors[tar, ] - exp(ln.mu + ln.var/2)) / sqrt((exp(ln.var) - 1) * exp(2*ln.mu + ln.var))
    } else if(errorInt == "gauss") {
      errors[tar, ] <- matrix(rnorm(length(tar) * n, mean = 0, sd = 1), nrow = length(tar), ncol = n)
    } else if (errorInt == "zero"){
      cst <- runif(n=length(tar), min=-3, max=3)
      errors[tar, ] <- matrix(cst, nrow = length(tar), ncol = n)
    }

    Y.k <- solve(diag(rep(1, v)) - B.k, errors)
    #max(abs(cov(t(Y.k))-sigma.k))
    data[[k]] <- Y.k
    Y <- cbind(Y, Y.k)
    Sigmas[[k]] <- sigma.k
  }

  return(list(Y = Y, data = data, B = B, Omega = Omega, B.true = B.true,
              Omega.true = Omega.true, targets = targets,
              target.length = target.length, Sigmas = Sigmas))
}


generateGraph <- function(v, n=NULL, l, d, b=0, scc_constraint=T, acy_constraint=F){
  # Setup B and Omega matrices
  B <- matrix(0, nrow = v, ncol = v)
  Omega <- diag(rep(1, v))
  
  # Generate cycle from 1 -> 2 -> ... K -> 1
  if(l > 1){
    for(i in 1:(l-1))
    {
      B[i+1, i] <- 1
    }
    B[1, l] <- 1
  }
  
  # fill in remaining edges
  for(j in 1:(v-1))
  {
    for(i in max(j+1, l+1):v)
    {
      if((B[i,j] != 1) & (B[j,i]!= 1)){
        U <- runif(1)
        if(U < d){
          B[i, j] <- 1
        } else if(U < (b + d))
        {
          Omega[i,j] <- Omega[j,i] <- 1
        }
      }
    }
  }
  
  
  # reorder the vertices
  reorder <- sample(v)
  B <- B[reorder, reorder]
  Omega <- Omega[reorder, reorder]
  
  return (list(B = B, Omega = Omega))
}

generateData <- function(B, Omega, targets=NULL, target.length=NULL, errorDist = "gauss", errorInt="gauss", det_constraint=0.01){
  # repeat until det(I-B.true) > det_constraint
  repeat {
    # Sample edge weights as given in the paper
    B.true <- matrix(runif(v^2,0.3,1), nrow = v) * B
    B.true <- matrix(rbinom(v^2,1,0.5)*2-1, nrow = v) * B.true
    
    if (det(diag(v)-B.true) > det_constraint)
      break
  }
  
  Omega.true <- diag(runif(v, 0.3, 1))
  
  
  data <- list()
  Y <- NULL
  Sigmas <- list()
  if (is.null(targets) || is.null(target.length)){
    targets = list(numeric(0))
    target.length = (n)
  }
  
  for (k in 1:length(target.length)){
    tar <- targets[[k]]
    n <- target.length[k]
    B.k <- B.true
    B.k[tar, ] <- 0
    Omega.k <- Omega.true
    Omega.k[tar, ] <- Omega.k[, tar] <- 0
    for (t in tar){
      for (i in 1:v){
        if (Omega[t, i] != 0 && !any(tar == i)){
          Omega.k[i, i] <- Omega.k[i, i] + sqrt(abs(Omega.true[t, i]*sqrt(Omega.true[i, i]*Omega.true[t, t])))
        }
      }
    }
    diag(Omega.k)[tar] <- diag(Omega.true)[tar]
    
    # Sample data from multivariate normal and make mean 0
    if(errorDist == "gamma")
    {
      errors <- matrix(rgamma(v * n, shape = 1, scale = 1), nrow = v, ncol = n)
    } else if(errorDist == "poisson")
    {
      pois.mu <- 3
      errors <- matrix(rpois(v * n, lambda = pois.mu), nrow = v, ncol = n)
      errors <- (errors - pois.mu)/sqrt(pois.mu)
    } else if(errorDist == "t"){
      t.df <- 5
      errors <- matrix(rt(v * n, df = t.df) * (t.df - 2)/t.df, nrow = v, ncol = n)
    } else if(errorDist == "lognormal"){
      ln.var <- 1
      ln.mu <- 0
      errors <- matrix(rlnorm(v * n, meanlog = ln.mu, sdlog = sqrt(ln.var)), nrow = v, ncol = n)
      errors <- (errors - exp(ln.mu + ln.var/2)) / sqrt((exp(ln.var) - 1) * exp(2*ln.mu + ln.var))
    } else if(errorDist == "gauss") {
      errors <- matrix(rnorm(v * n, mean = 0, sd = 1), nrow = v, ncol = n)
      #Y <- t(MASS::mvrnorm(n = n, mu = rep(0, v), Sigma = sigma))
    }
    
    #cov(t(errors))
    temp <- solve(diag(rep(1,v)) - B.k)
    sigma.k <- temp %*% Omega.k %*% t(temp)
    errors <- t(chol(Omega.k)) %*% errors
    
    # interventional errors distribution
    if(errorInt == "gamma")
    {
      errors <- matrix(rgamma(length(tar) * n, shape = 1, scale = 1), nrow = length(tar), ncol = n)
    } else if(errorInt == "poisson")
    {
      pois.mu <- 3
      errors[tar, ] <- matrix(rpois(length(tar) * n, lambda = pois.mu), nrow = length(tar), ncol = n)
      errors <- (errors - pois.mu)/sqrt(pois.mu)
    } else if(errorInt == "t"){
      t.df <- 5
      errors[tar, ] <- matrix(rt(length(tar) * n, df = t.df) * (t.df - 2)/t.df, nrow = length(tar), ncol = n)
    } else if(errorInt == "lognormal"){
      ln.var <- 1
      ln.mu <- 0
      errors[tar, ] <- matrix(rlnorm(length(tar) * n, meanlog = ln.mu, sdlog = sqrt(ln.var)), nrow = length(tar), ncol = n)
      errors[tar, ] <- (errors[tar, ] - exp(ln.mu + ln.var/2)) / sqrt((exp(ln.var) - 1) * exp(2*ln.mu + ln.var))
    } else if(errorInt == "gauss") {
      errors[tar, ] <- matrix(rnorm(length(tar) * n, mean = 0, sd = 1), nrow = length(tar), ncol = n)
    } else if (errorInt == "zero"){
      cst <- runif(n=length(tar), min=-3, max=3)
      errors[tar, ] <- matrix(cst, nrow = length(tar), ncol = n)
    }
    
    Y.k <- solve(diag(rep(1, v)) - B.k, errors)
    #max(abs(cov(t(Y.k))-sigma.k))
    data[[k]] <- Y.k
    Y <- cbind(Y, Y.k)
    Sigmas[[k]] <- sigma.k
  }
  
  return(list(Y = Y, data = data, B = B, Omega = Omega, B.true = B.true,
              Omega.true = Omega.true, targets = targets,
              target.length = target.length, Sigmas = Sigmas))
}

plotGraph <- function(
  g,
  tcltk=FALSE,
  ...
) ggm::plotGraph(g, layout=igraph::layout.circle, tcltk=tcltk, directed=T,...)


split_amat <- function (x) {
  res = list()
  B = x
  O = x
  B[B == 100] = 0
  O[O == 1] = 0
  O[O == 100] = 1
  diag(O) = 1
  res[[1]] = as.matrix(B)
  res[[2]] = as.matrix(O)
  return (res)
}

merge_amat <- function (B, O){
  if (any(dim(B) != dim(O)))
    stop("B and O must be of the same size!")
  graph = B + (O - diag(nrow(O))) * 100
  return (graph)
}