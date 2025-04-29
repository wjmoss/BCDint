ricf_R_int <- function(L = NULL, data, targets=NULL, target.length=NULL, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6),
                       maxiter=5000, out="none", maxkap = 1e13, B = NULL){
  
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
  S <- cov(data)
  
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
  tarjan <- function(u) {
    dft[u] <<- low[u] <<- count
    count <<- count + 1
    visited[u] <<- 1
    Stack <<- append(Stack, u)
    for (v in which(L[u, ] != 0)) {
      if (! visited[v]) {
        tarjan(v)
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
  
  # determine the nodes that only need 1 iteration
  # not in a cycle
  visited = rep(0, p)
  dft = rep(0, p)
  low = rep(0, p)
  count = 1
  Stack = list()
  for (u in 1:p) {
    if (! visited[u]) {
      tarjan(u)
    }
  }
  
  ## determine n1 (int out of C_k), n2 (int in C_k) for every node
  # determine pa and cc and intersection
  # check graph
  n1_i <- rep(0, p)
  n2_i <- rep(0, p)
  #invalid_i <- list()
  posind <- cumsum(c(0,target.length))
  cov_list <- list()
  cc_i <- list()
  pa_i <- list()
  j_i <- list()
  for (i in 1:p){
    #n1_i[i] <- nrow(data) %/% 2
    #n2_i[i] <- nrow(data) - n1_i[i]
    #rows_i[[i]] <- c(1:nrow(data))
    cc_i[[i]] <- which(low == low[i])
    pa_i[[i]] <- which(L[, i] != 0)
    #tmp <- intersect(cc_i[[i]], pa_i[[i]])
    pa_cc <- which(L[, i] != 0 & low == low[i])
    if (length(pa_cc) > 1){
      stop("The graph is invalid for the algorithm!")
    }
    else {
      if (length(pa_cc) == 0)
        j_i[[i]] <- numeric(0)
      else
        j_i[[i]] <- pa_cc
    }
    #j_i[[i]] <- intersect(cc_i[[i]], pa_i[[i]])
    invalid_tar <- which(sapply(targets, function(x) any(x == i)))
    if (length(invalid_tar) == length(targets))
      stop(paste("No observational data for node", i, "!", sep=" "))
    ind <- c()
    for (tar in invalid_tar){
      ind <- c(ind, (posind[tar]+1): posind[tar+1])
    }
    #invalid_rows <- which(t$target.index %in% invalid_targets)
    #rows_i[[i]] <- setdiff(c(1:nrow(data)), invalid_rows)
    if (length(ind) == 0){
      cov_list[[i]] <- cov(data)
    }
    else{
      cov_list[[i]] <- cov(data[-ind, ])
    }
    n1_ind <- which(sapply(targets, function(x) length(intersect(cc_i[[i]], x))==0))
    n1_i[i] <- sum(target.length[n1_ind])
    n2_i[i] <- nrow(data) - n1_i[i] - sum(target.length[invalid_tar])
    #n2[i] <- nrow(data) - length(invalid_rows) - n1[i]
  }
  
  
  # ricf_dg works for n1=0 or n2=0
  # bcd_int version, including root selection
  
  Lcur <- Linit; Ocur <- Oinit
  iter <- 1
  repeat {
    
    if(p>1){
      for (i in 1:p) {
        if (iter > 1) {
          if (dft[i] == low[i] && visited[i] == 1) next
        }
        #cc <- which(low == low[i])
        #pa <- which(L[, i] != 0)
        n1 <- n1_i[i]
        n2 <- n2_i[i]
        S <- cov_list[[i]]
        cc <- cc_i[[i]]
        pa <- pa_i[[i]]
        n.pa <- length(pa)
        j <- j_i[[i]]
        #IB <- diag(p) - Lcur
        #Elessi <- IB[, -i]
        
        if (max(Ocur[-i]) / min(Ocur[-i]) > maxkap) {
          stop(paste("The condition number of Ocur[-i, -i] is too large for node", i))
        }
        # elementwise product == row product
        #Zlessi <- t(t(Elessi) * Ocur[-i])
        #Zhelp <- S %*% Zlessi
        # The following line gets the indices of Zlessi corresponding to spouses
        #zsp <- c(sp[sp < i], sp[sp > i] - 1)
        
        # n1 or n2 = 0
        # to do:simpler version?
        if (n1 == 0 || n2 == 0){
          #S = cov(data[-delrows, , drop=FALSE])
          if (n.pa > 0) {
            ## DETERMINE a AND a0
            ## a: coefficient of B[i,pa(i)]
            ## a0: const independent to B[i,pa(i)]
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
            #pa.pos <- pa[ind.pos]
            ## general case
            if (length(ind.pos) > 0) {
              ## can be simplied if no siblings
              M <- S[pa, pa]
              m <- S[pa, i]
              alpha <- solve(S[pa, pa], S[pa, i])
              y0 <- S[i, i] - S[i, pa] %*% alpha
              coef <- alpha + as.numeric(y0 / (a0 + t(a) %*% alpha)) * solve(S[pa, pa], a)
              Lcur[pa, i] <- coef
              ## n.sp > 0 case
            }
            
            ## a = 0, reduced to Linear Regression
            else {
              ## Perform linear regression for Betas and (possibly) Omegas
              M <- S[pa, pa]
              m <- S[pa, i]
              coef <- solve(M, m)
              y0 <- S[i, i] - t(coef) %*% m
              Lcur[pa, i] <- coef
              
            }
            
            
            ## Find the variance omega_{ii}
            RSS_avg <- S[i, i] - 2 * t(as.matrix(coef)) %*% m + t(as.matrix(coef)) %*% M %*% as.matrix(coef)
            if (max(Ocur[-i]) / min(Ocur[-i]) > maxkap) {
              stop(paste("The Condition number of Ocur[-i, -i] is too large for node", i))
            }
            
            ### print ###
            
            Ocur[i] <- RSS_avg
          }
          else{
            Ocur[i] <- S[i, i]
          }
        }
        # n1, n2 != 0
        else{
          #compute auxiliary values lambda,a,b,c
          tmp <- 1 - det((diag(p) - Lcur)[cc, cc, drop=FALSE])
          Lcur[j, i] <- Lcur[j, i] + 1
          Clessi <- 1 - det((diag(p) - Lcur)[cc, cc, drop=FALSE]) - tmp
          Lcur[j, i] <- Lcur[j, i] - 1
          pa_adj <- c(pa[pa < j], pa[pa > j])
          if (length(pa_adj) > 0){
            mat_tmp <- solve(S[pa_adj, pa_adj])
            a <- S[j, j] - S[j, pa_adj] %*% mat_tmp %*% S[pa_adj, j]
            b <- -(S[j, i] - S[j, pa_adj] %*% mat_tmp %*% S[pa_adj, i])
            c <- S[i, i] - S[i, pa_adj] %*% mat_tmp %*% S[pa_adj, i]
          }
          else{
            a <- S[j, j]
            b <- -S[j, i]
            c <- S[i, i]
          }
          
          # change to ratio of n1 and n2?
          #delta <- (b*Clessi*(n1-n2) + a*(n1+n2))^2 + 4*a*Clessi*n2*(c*Clessi*n1+b*(n1+n2))
          ratio <- n1 / n2
          delta <- (b*Clessi*(ratio-1) + a*(ratio+1))^2 + 4*a*Clessi*(c*Clessi*ratio+b*(ratio+1))
          if (tmp == 1)
            stop("Singular value of edge weights!")
          # compare two possible solutions
          coef_1 <- coef_2 <- rep(0, nrow(L))
          coef_1[j] <- ( a*(ratio+1) + b*Clessi*(ratio-1) - sqrt(delta) ) / (2*a*Clessi)
          coef_2[j] <- ( a*(ratio+1) + b*Clessi*(ratio-1) + sqrt(delta) ) / (2*a*Clessi)
          
          #if (Clessi > 0)
          # cycle -i prod > 0, choose x<1/l
          #coef_ji <- ( a*(ratio+1) + b*Clessi*(ratio-1) - sqrt(delta) ) / (2*a*Clessi)
          #else
          # cycle -i prod < 0, choose x>1/l
          #coef_ji <- ( a*(ratio+1) + b*Clessi*(ratio-1) + sqrt(delta) ) / (2*a*Clessi)
          # Clessi = 0?
          
          
          # choose different branches, Lcur[j,i]<1/Clessi ==> left branch?
          # avoid likelihood comparison
          #if (tmp == 1)
          #stop("Singular value of edge weights!")
          #if (Clessi > 0)
          # cycle -i prod > 0, choose x<1/l
          #coef_ji <- ( a*(n1+n2) + b*Clessi*(n1-n2) - sqrt(delta) ) / (2*a*Clessi*n2)
          #else
          # cycle -i prod < 0, choose x>1/l
          #coef_ji <- ( a*(n1+n2) + b*Clessi*(n1-n2) + sqrt(delta) ) / (2*a*Clessi*n2)
          # Clessi = 0?
          
          # compare two log-likelihood values
          # difference between solve(A) * b and solve(A,b)?
          if (length(pa_adj) != 0){
            coef_1[pa_adj] <- mat_tmp %*% (S[pa_adj, i] - S[pa_adj, j] * coef_1[j])
          }
          omega_i_1 <- S[i, i] - 2 * t(coef_1[pa]) %*% S[pa, i] + t(coef_1[pa]) %*% S[pa, pa] %*% coef_1[pa]
          loglike_1 <- ratio * log((1-coef_1[j]*Clessi)^2) - (ratio+1)*log(omega_i_1)
          
          if (length(pa_adj) != 0){
            coef_2[pa_adj] <- mat_tmp %*% (S[pa_adj, i] - S[pa_adj, j] * coef_2[j])
          }
          omega_i_2 <- S[i, i] - 2 * t(coef_2[pa]) %*% S[pa, i] + t(coef_2[pa]) %*% S[pa, pa] %*% coef_2[pa]
          loglike_2 <- ratio * log((1-coef_2[j]*Clessi)^2) - (ratio+1)*log(omega_i_2)
          
          if (is.na(loglike_1) || is.na(loglike_2))
            stop("NA value of log-likelihood!")
          if (loglike_1 > loglike_2){
            Lcur[pa, i] <- coef_1[pa]
            Ocur[i] <- omega_i_1
          } else{
            Lcur[pa, i] <- coef_2[pa]
            Ocur[i] <- omega_i_2
          }
          
          
          ## Find the variance omega_{ii}
          #if (length(pa_adj) != 0)
          #Lcur[pa_adj, i] <- mat_tmp %*% (S[pa_adj, i] - S[pa_adj, j] * coef_ji)
          #Lcur[j, i] <- coef_ji
          #Ocur[i] <- S[i, i] - 2 * t(Lcur[pa, i]) %*% S[pa, i] + t(Lcur[pa, i]) %*% S[pa, pa] %*% Lcur[pa, i]
          if (max(Ocur[-i]) / min(Ocur[-i]) > maxkap) {
            stop(paste("The Condition number of Ocur[-i, -i] is too large for node", i))
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
      
      #rows <- rows_i[[1]]
      S = cov(data[rows, , drop=FALSE])
      sigcur = as.matrix(S[p, p])
      bhat = as.matrix(0)
      Ocur = S[p, p]
      Lcur = as.matrix(0)
      iter = 1
      break
    }
    
  }
  return(list(Sigmahat = sigcur, Bhat = bhat, Omegahat = Ocur, Lambdahat = Lcur, iterations = iter, converged = (iter < maxiter)))
}
