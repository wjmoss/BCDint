library(purrr)
library(dplyr)
library(tidyverse)
library(MASS)
library(ggplot2)
library(SEMgraph)
source("ricf_int.R")
source("ricf_dg.R")

## get the relevant data and target.length
names(sachs$rawdata)
idx <- c(1, 3:6)
target.length <- sapply(idx, \(i) nrow(sachs$rawdata[[i]]))
sum(target.length)

## take log-transformation, and merge the relevant data into a large matrix of sum(N_k) * 11
combined <- bind_rows(lapply(sachs$rawdata[idx], function(df) {
  df %>%
    mutate(across(where(is.numeric), log)) %>%          # log-transform numeric columns
    mutate(across(where(is.numeric), ~ .x - mean(.x))) # center numeric columns
}))

## check dim and visualization of combined data frame
dim(combined)

data_long <- combined %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# one ggplot with facets
ggplot(data_long, aes(x = value)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of each variable", x = "Value", y = "Count")




## define the graph and targets
L <-  matrix(0, 11, 11)
rownames(L) <- colnames(L) <- names(combined)

# acyclic graph
L[1, 2] <- 1
L[2, 6] <- 1
L[3, c(4,9)] <- 1
L[4, 9] <- 1
L[5, c(3,4,7)] <- 1
L[6, 7] <- 1
L[8, c(1,2,6,7,10,11)] <- 1
L[9, c(1,2,8,10,11)] <- 1
sum(L)

## data1: null ()
## data3: akt (7)
## data4: pkc (9)
## data5: pip2 (4)
## data6: mek (2)
targets <- list(numeric(0), c(7), c(9), c(4), c(2))


## run BCDint
res_a <- ricf_int(L = L,
                 data = as.matrix(combined), 
                 targets = targets,
                 target.length = target.length)

# cyclic graph
L[5, 4] <- 0
L[4, 5] <- 1
res_c <- ricf_int_(L = L,
                 data = as.matrix(combined), 
                 targets = targets,
                 target.length = target.length)

## compute and compare the likelihoods
compute_llh <- function(res, data, targets, target.length){
  val <- 0
  
  v <- ncol(data)
  ind <- c(0, cumsum(target.length))
  for (k in seq_along(targets)){
    target <- targets[[k]]
    L <- res$Lambdahat
    L[, target] <- 0
    O <- res$Omegahat
    
    data_tmp <- data[(ind[k]+1):ind[k+1], ]
    for (i in target){
      O[i] <- var(data_tmp[, i])
    }
    
    Sigma <- cov(data_tmp)
    K <- (diag(v) - res$Lambdahat) %*% (1 / O * t(diag(v)-L))
    tr <- sum( c(K) * c(t(Sigma)) )
    val <- val + (-sum(log(O)) + log(det(diag(v)-L) ^ 2) - tr) * target.length[k]
  }
  
  return (val)
}

compute_llh(res_a, combined, targets, target.length)
compute_llh(res_c, combined, targets, target.length)


## use the directed part in Drton (2019)
L <-  matrix(0, 11, 11)
rownames(L) <- colnames(L) <- names(combined)

# acyclic graph
L[1, 2] <- 1
L[2, 6] <- 1
L[3, c(4,9)] <- 1
L[4, 9] <- 1
L[5, c(3,4,7)] <- 1
L[8, c(1,6,7,10,11)] <- 1
L[9, c(1,10,11)] <- 1
res_a <- ricf_int(L = L,
                  data = as.matrix(combined), 
                  targets = targets,
                  target.length = target.length)

# cyclic graph
L[5, 4] <- 0
L[4, 5] <- 1
res_c <- ricf_int_(L = L,
                   data = as.matrix(combined), 
                   targets = targets,
                   target.length = target.length)

compute_llh(res_a, combined, targets, target.length)
compute_llh(res_c, combined, targets, target.length)

