## RICF for interventional data, for arbitrary digraph and special interventions

# INPUTS ---
# L, O: the 0-1 valued adjacency matrices indicating the structure of the graph G
# data: n-by-p matrix of observational data and interventional data
# targets: list, unique intervention targets
# target.length: vector, the length of every interventional data set, in the order of targets
# Linit; Oinit: the initial values of the edge parameters for running the algorithm
# sigconv: Boolean (True = look for conv. in Sigma, False = look for conv. in L/O )
# tol: gives the max average entrywise absolute deviation in Sigma permited for convergence
# maxiter: gives the max number of iterations before divergence is accepted
# out: String: options: None/False, Final, All/True
# maxkap: maximum condition number accepted before error thrown
# B: (optional instead of L -- here B = t(L))
# covariance: unbiased uses cov(); ml uses Y'Y/n for already centered/zero-mean data
#
# OUTPUTS ---
# Sigmahat: the fitted value for Sigma resulting from the algorithm
# Bhat: the fitted value for the directed edge matrix (note: equal to I-t(Lambdahat))
# Omegahat, Lambdahat: the MLE for the parameter values resulting from the algorithm
# iterations: the number of iterations run by the algorithm before convergence or divergence accepted
# converged: TRUE or FALSE - based on whether or not the algorithm converged before maxiter reached
ricf_int_ <- function(L = NULL, data, targets=NULL, target.length=NULL, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6),
                      maxiter=5000, out="none", maxkap = 1e13, B = NULL,
                      covariance=c("unbiased", "ml")){

  covariance <- match.arg(covariance)
  sample_moment <- function(Y) {
    if (covariance == "ml") crossprod(Y) / nrow(Y) else cov(Y)
  }

  ## parameter validations
  # matrices
  if (is.null(L)) {
    if (!is.matrix(B))
      stop("B must be a matrix!")
    if (nrow(B) != ncol(B))
      stop("B must be a square matrix!")
    L = t(B)
  }
  if (!is.matrix(L))
    stop("L must be a matrix!")
  if (!is.matrix(data))
    stop("data must be a matrix!")
  if (ncol(data) != ncol(L))
    stop("ncol of data must equal to the number of variables!")
  if (nrow(data) < nrow(L))
    stop("nrow of data must be at least the number of variables!")
  if (nrow(L) != ncol(L))
    stop("L must be a square matrix!")
  p <- nrow(L)
  O <- rep(1, p)
  S <- sample_moment(data)

  #targets
  if (is.null(targets) || is.null(target.length)){
    #all observation data
    targets = list(numeric(0))
    target.length = nrow(data)
  }
  else{
    if (!all(sapply(targets, is.vector)) || !all(sapply(targets, is.numeric))){
      stop("Every target must be a numeric vector!")
    }
  }
  if (any(sapply(targets, function(x) (length(x) != 0 && (1 > min(x) || max(x) > p))) ))
    stop("Invalid target!")
  if (!is.null(target.length) && (sum(target.length) != nrow(data)))
    stop("Wrong data counts!")
  # preprocessing of intervention targets
  #targets = lapply(targets, FUN = unique)

  #if (length(target) == 0 || (1 <= min(target) && max(target) <= p))
    #stop("Invalid target!")

  # Initialize the directed edge parameters via OLS
  initL <- function(L, S) {
    Linit <- L
    parents <- apply(L, 2, function(x) which(x > 0))
    len <- apply(L, 2, sum)
    for (j in 1:p) {
      if (len[j] > 0) {
        for (k in 1:len[j]) {
          Linit[parents[[j]][k], j] <- S[parents[[j]][k], j] / S[parents[[j]][k], parents[[j]][k]]
          # 2-cycle adjustment
          if (L[j, parents[[j]][k]] == 1)
            Linit[parents[[j]][k], j] <- Linit[parents[[j]][k], j] / 2
        }
      }
    }
    return(Linit)
  }

  # Initialize the variances by diag(S)
  initO <- function(S) {
    ## O2 <- diag(diag(S))  ## for BAPs
    ## R <- diag(p)
    ## if(p>1)
    ##   for(i in 1:(p-1)) { for(j in (i+1):p) { R[j, i] <- R[i, j] <- rnorm(1) }}
    ## O2 <- R*O
    ## diag(O2) <- rep(0, p)
    ## diag(O2) <- rowSums(abs(O2)) + abs(rnorm(p))
    return(diag(S))
  }

  if (is.null(Linit)) { Linit <- initL(L, S) }
  if (is.null(Oinit)) { Oinit <- initO(S) }
  if (!is.matrix(Linit))
    stop("Linit need to be matrix!")
  if (!is.matrix(Oinit) && !is.vector(Oinit)){
    stop("Oinit need to be matrix or vector!")
  }
  if (is.matrix(Oinit)){
    Oinit <- diag(Oinit)
  }
  if (!all(Oinit > 0)){
    stop("All values in Oinit must be positive!")
  }
  if (nrow(Linit) != ncol(Linit))
    stop("Linit must be square matrix!")
  if (length(unique(c(nrow(Linit), length(Oinit), nrow(L), ncol(S)))) > 1)
    stop("One of the input matrices has the wrong dimension!")
  if (maxiter <= 0 || maxiter %% 1 != 0)
    stop("A positive integer is needed for the max number of iterations!")
  if (tol <= 0)
    stop("A positive tolerance is needed for convergence to be possible!")
  if (!is.logical(sigconv))
    stop("sigconv needs to take on a logical value!")
  out = tolower(out)
  if (!is.character(out))
    stop("Output needs to be a string: none/false, final, or all/true!")
  if (!(out == "true" || out == "all" || out == "final" || out == "false" || out == "none")) {
    stop("Output variable needs to be: none/false, final, or all/true!")
  }

  Det <- function(Lcur) {
    return(det(diag(nrow(Lcur)) - Lcur))
  }

  ## Tarjan's algorithm, to find all strongly connected components
  tarjan <- function(u, L) {
    dft[u] <<- low[u] <<- count
    count <<- count + 1
    visited[u] <<- 1
    Stack <<- append(Stack, u)
    for (v in which(L[u, ] != 0)) {
      if (! visited[v]) {
        tarjan(v, L)
        low[u] <<- min(low[u], low[v])
      } else if (v %in% Stack) {
        low[u] <<- min(low[u], dft[v])
      }
    }
    if (dft[u] == low[u]) {
      while(length(Stack) > 0) {
        if (Stack[[length(Stack)]] != u) {
          Stack[[length(Stack)]] <<- NULL
          visited[u] <<- visited[u] + 1
        } else {
          Stack[[length(Stack)]] <<- NULL
          break
        }
      }
    }
  }

  visited = rep(0, p)
  dft = rep(0, p)
  low = rep(0, p)
  count = 1
  Stack = list()
  for (u in 1:p) {
    if (! visited[u]) {
      tarjan(u, L)
    }
  }

  # determine the nodes that only need 1 iteration: single node strongly connected component {i}
  one_iter <- (dft == low & visited == 1)
  # equivalent
  #G <- graph_from_adjacency_matrix(L)
  #cps <- components(G, mode="strong")
  #one_iter <- cps$csize[cps$membership[1:v]] == 1

  ## determine n1 (scc = scc_unique), n2 (scc = {i}) for every node
  ## check interventions
  n1_i <- rep(0, p)
  n2_i <- rep(0, p)
  posind <- cumsum(c(0,target.length))
  cov_list <- list()
  for (i in 1:p){
    # find valid intervention targets for node i
    invalid_tar <- which(sapply(targets, function(x) any(x == i)))
    if (length(invalid_tar) == length(targets))
      stop(paste("No observational data for node", i, "!", sep=" "))
    ind <- integer()
    for (tar in invalid_tar){
      # invalid indices for node i, to remove
      ind <- c(ind, (posind[tar]+1): posind[tar+1])
    }

    # compute the covariance matrix w.r.t. valid data rows
    if (length(ind) == 0){
      cov_list[[i]] <- sample_moment(data)
    } else{
      cov_list[[i]] <- sample_moment(data[-ind, , drop=FALSE])
    }

    # check scc in intervened graphs and compute n1, n2
    cp_unique <- NULL
    valid_tar <- setdiff(1:length(targets), invalid_tar)
    for (i_tar in valid_tar){
      # classify the intervention to two groups
      tar <- targets[[i_tar]]
      Ltar <- L
      Ltar[tar, ] <- 0
      
      # check strongly connected components
      #too long time!
      #Gtar <- graph_from_adjacency_matrix(Ltar)
      #cp_nodes <- intersect(subcomponent(Gtar, i, mode="in"), subcomponent(Gtar, i, mode="out"))
      visited = rep(0, p)
      dft = rep(0, p)
      low = rep(0, p)
      count = 1
      Stack = list()
      tarjan(i, Ltar)
      cp_nodes <- (1:p)[low == low[i]]
      
      if (length(cp_nodes) == 1)
        n2_i[i] <- n2_i[i] + target.length[i_tar]
      else{
        cp <- Ltar[cp_nodes, cp_nodes]
        if (is.null(cp_unique)){
          cp_unique <- cp
          n1_i[i] <- n1_i[i] + target.length[i_tar]
        }
        else{
          if (length(cp) == length(cp_unique) && all(cp == cp_unique))
            n1_i[i] <- n1_i[i] + target.length[i_tar]
          else
            stop("Two different non-trivial strongly connected componets containing i, BCD algorithm has degree >= 3!")
        }
      }
    }

  }


  Lcur <- Linit; Ocur <- Oinit
  iter <- 1
  repeat {
    #print(iter)
    #print(Lcur)
    #print(Ocur)
    if(p>1){
      for (i in 1:p) {
        if (iter > 1) {
          if (one_iter[i] == 1) next
        }
        n1 <- n1_i[i]
        n2 <- n2_i[i]
        S <- cov_list[[i]]
        pa <- which(L[, i] != 0)
        n.pa <- length(pa)

        if (max(Ocur[-i]) / min(Ocur[-i]) > maxkap) {
          stop(paste("The condition number of Ocur[-i, -i] is too large for node", i))
        }

        if (n.pa == 0){
          Ocur[i] <- S[i, i]
        }
        else{
          a <- rep(0, n.pa)
          for (k in 1:n.pa) {
            temp <- Lcur
            temp[pa[k], i] <- 2
            det2 <- Det(temp)
            temp[pa[k], i] <- 1
            det1 <- Det(temp)
            a[k] <- det2 - det1
          }
          temp[pa, i] <- 0
          a0 <- Det(temp)

          ind.pos <- which(a != 0)
          alpha <- solve(S[pa, pa], S[pa, i])
          y0 <- S[i, i] - S[i, pa] %*% alpha # y0^2 in text

          if (n1 == 0 || length(ind.pos) == 0){
            # reduced to linear regression
            #M <- S[pa, pa]
            #m <- S[pa, i]
            Lcur[pa, i] <- alpha
            Ocur[i] <- S[i, i] - 2 * t(as.matrix(alpha)) %*% S[pa, i] + t(as.matrix(alpha)) %*% S[pa, pa] %*% as.matrix(alpha)
          }
          else{
            if (n2 == 0){
              # equi to obs case, deg=1
              #M <- S[pa, pa]
              #m <- S[pa, i]
              coef <- alpha + as.numeric(y0 / (a0 + t(a) %*% alpha)) * solve(S[pa, pa], a)
              Lcur[pa, i] <- coef
              Ocur[i] <- S[i, i] - 2 * t(as.matrix(coef)) %*% S[pa, i] + t(as.matrix(coef)) %*% S[pa, pa] %*% as.matrix(coef)
            }
            else{
              # our algorithm, deg=2
              # compute auxiliary quantities
              ratio <- n1 / n2
              aa <- sum(a * solve(S[pa, pa], a))
              bb <- (sum(a * alpha) + a0) * (ratio + 1)
              cc <- - ratio * y0

              # two possible solutions and log_likelihood values
              t1 <- (- bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
              t2 <- (- bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
              coef_1 <- alpha + t1 * solve(S[pa, pa], a)
              coef_2 <- alpha + t2 * solve(S[pa, pa], a)
              omega_1 <- S[i, i] - 2 * t(as.matrix(coef_1)) %*% S[pa, i] + t(as.matrix(coef_1)) %*% S[pa, pa] %*% as.matrix(coef_1)
              omega_2 <- S[i, i] - 2 * t(as.matrix(coef_2)) %*% S[pa, i] + t(as.matrix(coef_2)) %*% S[pa, pa] %*% as.matrix(coef_2)
              loglike_1 <- ratio * log((a0+sum(coef_1*a))^2) - (ratio+1)*log(omega_1)
              loglike_2 <- ratio * log((a0+sum(coef_2*a))^2) - (ratio+1)*log(omega_2)

              # compare the log-likelihood values
              if (is.na(loglike_1) || is.na(loglike_2))
                stop("NA value of log-likelihood!")
              if (loglike_1 > loglike_2){
                Lcur[pa, i] <- coef_1
                Ocur[i] <- omega_1
              } else{
                Lcur[pa, i] <- coef_2
                Ocur[i] <- omega_2
              }
            }
          }

        }

      }


      #sigcur <- t(solve(diag(p) - Lcur)) %*% diag(Ocur) %*% solve(diag(p) - Lcur)
      sigcur <- t(solve(diag(p) - Lcur)) %*% (Ocur * solve(diag(p) - Lcur))
      bhat <- diag(p)-t(Lcur)
      if (iter == 1) {
        if (out == "true" || out == "all"){cat(iter, "\n")}
        if (maxiter == 1) {
          if (out == "true" || out == "all" || out == "final") {
            cat("Sigmahat \n")
            print(sigcur)
            cat("\nBhat \n")
            print(bhat)
            cat("\nOmegahat \n")
            print(Ocur)
            cat("\nLambdahat \n")
            print(Lcur)
            cat("\niterations \n")
            print(iter)
          }
          break
        }
      }
      else if (iter > 1) {
        dsig <- mean(abs(sigcur - sigpast))
        dLO <- sum(abs(Lcur - Lpast) + abs(Ocur - Opast)) / (sum(L) + sum(O))
        if (out == "true" || out == "all") {
          dsig6 <- format(round(dsig, 6), nsmall = 6)
          dLO6 <- format(round(dLO, 6), nsmall = 6)
          cat(iter, " Avg Change in L & O: ", dLO6, "| Avg Change in Sigma: ", dsig6, "\n")
        }

        if ((sigconv && dsig < tol) || (!sigconv && dLO < tol) || iter >= maxiter) {
          if (out == "true" || out == "all" || out == "final") {
            cat("Sigmahat \n")
            print(sigcur)
            cat("\nBhat \n")
            print(bhat)
            cat("\nOmegahat \n")
            print(Ocur)
            cat("\nLambdahat \n")
            print(Lcur)
            cat("\niterations \n")
            print(iter)
          }
          break
        }
      }
      sigpast <- sigcur; Lpast <- Lcur; Opast <- Ocur
      iter <- iter + 1
    }

    else {

      S = cov_list[[1]]
      if (!is.finite(S[1, 1]) || S[1, 1] <= 0)
        stop("Single-variable fitting requires positive variance in valid data")
      sigcur = as.matrix(S[p, p])
      bhat = as.matrix(1)
      Ocur = S[p, p]
      Lcur = as.matrix(0)
      iter = 1
      break
    }

  }
  return(list(Sigmahat = sigcur, Bhat = bhat, Omegahat = Ocur, Lambdahat = Lcur, iterations = iter, converged = (p == 1 || iter < maxiter)))
}


## compute the log-likelihood of a covariance matrix S
llh <- function(obj, S)
{
  K = (diag(nrow(S))-obj$Lambdahat) %*% (1 / obj$Omegahat * t(diag(nrow(S))-obj$Lambdahat))
  tr <- sum( c(K) * c(S) )
  llh <- - sum(log(obj$Omegahat)) + log(det(diag(nrow(S))-obj$Lambdahat) ^ 2) - tr
  return (llh)
}


# Environment-weighted score; constants and the 1/2 are omitted.
# ml uses Y'Y/n for zero-mean data; unbiased retains legacy cov() normalization.
llh_int <- function(obj, data, targets=NULL, target.length=NULL,
                    covariance=c("unbiased", "ml")) {
  covariance <- match.arg(covariance)
  data <- as.matrix(data)
  p <- ncol(data)
  if (is.null(targets) || is.null(target.length)) {
    targets <- list(numeric(0))
    target.length <- nrow(data)
  }
  if (length(targets) != length(target.length) ||
      any(!is.finite(target.length) | target.length < 2 | target.length %% 1 != 0) ||
      sum(target.length) != nrow(data))
    stop("Each environment must have at least two rows and matching counts")
  offsets <- c(0, cumsum(target.length))
  value <- 0
  for (k in seq_along(targets)) {
    target <- targets[[k]]
    L <- obj$Lambdahat
    L[, target] <- 0
    O <- obj$Omegahat
    Y <- data[seq.int(offsets[k]+1, offsets[k+1]), , drop=FALSE]
    S <- if (covariance == "ml") crossprod(Y) / nrow(Y) else cov(Y)
    O[target] <- diag(S)[target]
    if (any(!is.finite(O) | O <= 0)) return(-Inf)
    A <- diag(p) - L
    K <- A %*% (1 / O * t(A))
    value <- value + nrow(Y) * (-sum(log(O)) +
      2 * as.numeric(determinant(A, logarithm=TRUE)$modulus) - sum(K * t(S)))
  }
  value
}

# Compare starts using all environments, and forward the solver controls.
ricf_int <- function(L=NULL, data, targets=NULL, target.length=NULL,
                     Linit=NULL, Oinit=NULL, sigconv=TRUE, tol=1e-6,
                     maxiter=5000, out="none", maxkap=1e13, B=NULL, restarts=1,
                     covariance=c("unbiased", "ml")) {
  covariance <- match.arg(covariance)
  if (is.null(L)) {
    if (!is.matrix(B)) stop("Supply a matrix L or B")
    L <- t(B)
  }
  v <- nrow(L)
  if (length(restarts) != 1 || !is.finite(restarts) || restarts < 0 ||
      restarts %% 1 != 0) stop("restarts must be a nonnegative integer")
  starts <- lapply(seq_len(restarts), function(i)
    list(Linit=matrix(runif(v*v,-1,1),v,v)*L, Oinit=runif(v,0,1)))
  starts <- c(starts, list(list(Linit=L*0, Oinit=diag(v)),
                          list(Linit=Linit, Oinit=Oinit)))
  fits <- lapply(starts, function(start) ricf_int_(
    L=L, data=data, targets=targets, target.length=target.length,
    Linit=start$Linit, Oinit=start$Oinit, sigconv=sigconv, tol=tol,
    maxiter=maxiter, out=out, maxkap=maxkap, covariance=covariance))
  scores <- vapply(fits, llh_int, numeric(1), data=data,
                   targets=targets, target.length=target.length, covariance=covariance)
  if (!any(is.finite(scores))) stop("No finite joint score among initializations")
  selected <- which.max(scores)
  result <- fits[[selected]]
  result$restart_scores <- scores
  result$selected_restart <- selected
  result$joint_score <- scores[selected]
  result$covariance <- covariance
  result
}
