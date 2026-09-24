# function
f <- function(x, n1 = 1, n2 = 1,
              a, b, c, d, lambda) {
  
  den1 <- x + d / lambda
  den2 <- a * x^2 + 2 * b * x + c
  
  y <- n1 / den1 -
    (n1 + n2) * (a * x + b) / den2
  
  # remove points very close to singularities
  y[abs(den1) < 1e-3] <- NA
  y[abs(den2) < 1e-3] <- NA
  
  y
}


# find zeros
get_zeros <- function(n1 = 1, n2 = 1,
                      a, b, c, d, lambda,
                      tol = 1e-8) {
  
  # coefficients of A*x^2 + B*x + C = 0
  A <- n1 * a - (n1 + n2) * a
  
  B <- 2 * n1 * b -
    (n1 + n2) * (a * d / lambda + b)
  
  C <- n1 * c -
    (n1 + n2) * b * d / lambda
  
  # candidate roots
  roots <- numeric(0)
  
  if (abs(A) < tol) {
    
    if (abs(B) >= tol) {
      roots <- -C / B
    }
    
  } else {
    
    delta <- B^2 - 4 * A * C
    
    if (delta >= -tol) {
      
      # numerical protection
      delta <- max(delta, 0)
      
      roots <- c(
        (-B - sqrt(delta)) / (2 * A),
        (-B + sqrt(delta)) / (2 * A)
      )
    }
  }
  
  if (length(roots) == 0)
    return(numeric(0))
  
  # ------------------------------------------------
  # IMPORTANT:
  # remove roots outside the domain of original f(x)
  # ------------------------------------------------
  
  den1 <- roots + d / lambda
  den2 <- a * roots^2 + 2 * b * roots + c
  
  valid <- abs(den1) > tol &
    abs(den2) > tol
  
  roots <- roots[valid]
  
  # remove duplicated roots
  unique(round(roots, 10))
}

## to plot
# x range
x <- seq(-3, 7, length.out = 10000)


# four cases
pars <- list(
  
  list(
    type = "Type (i)",
    a = 0, b = 0, c = 2,
    d = -1, lambda = 1
  ),
  
  list(
    type = "Type (ii)",
    a = 1, b = 1, c = 2,
    d = -1, lambda = 1
  ),
  
  list(
    type = "Type (iii)",
    a = 1, b = 1, c = 1,
    d = -1, lambda = 1
  ),
  
  list(
    type = "Type (iv)",
    a = 1, b = 1, c = 1,
    d = -1, lambda = -1
  )
)


par(mfrow = c(2, 2),
    mar = c(4, 4, 4, 1))


for (p in pars) {
  
  y <- f(
    x,
    n1 = 1,
    n2 = 1,
    a = p$a,
    b = p$b,
    c = p$c,
    d = p$d,
    lambda = p$lambda
  )
  
  # Prevent huge values from messing up the graph scale
  y[abs(y) > 10] <- NA
  
  plot(
    x, y,
    type = "l",
    xlim = c(-3, 7),
    ylim = c(-4, 4),
    xlab = "x",
    ylab = "f(x)",
    main = p$type
  )
  
  # horizontal line y = 0
  abline(
    h = 0,
    col = "red",
    lty = 2
  )
  
  # singularity from x + d/lambda = 0
  pole1 <- -p$d / p$lambda
  
  abline(
    v = pole1,
    col = "blue",
    lty = 2
  )
  
  # roots of ax^2 + 2bx + c
  if (p$a != 0) {
    
    delta <- (2 * p$b)^2 -
      4 * p$a * p$c
    
    if (delta >= 0) {
      
      roots <- c(
        (-2*p$b + sqrt(delta)) / (2*p$a),
        (-2*p$b - sqrt(delta)) / (2*p$a)
      )
      
      for (r in unique(round(roots, 10))) {
        
        # 避免和已有 pole 重复画线
        if (abs(r - pole1) > 1e-8) {
          abline(
            v = r,
            col = "blue",
            lty = 2
          )
        }
      }
    }
  }
  
  # zeros
  zeros <- get_zeros(
    n1 = 1,
    n2 = 1,
    a = p$a,
    b = p$b,
    c = p$c,
    d = p$d,
    lambda = p$lambda
  )
  
  if (length(zeros) > 0) {
    
    # only show zeros in xlim
    zeros <- zeros[
      zeros >= -3 &
        zeros <= 7
    ]
    
    points(
      zeros,
      rep(0, length(zeros)),
      col = "red",
      pch = 19
    )
  }
  
  # parameter subtitle
  mtext(
    sprintf(
      "a = %g, b = %g, c = %g, d = %g, lambda = %g",
      p$a, p$b, p$c, p$d, p$lambda
    ),
    side = 3,
    line = 0.3,
    cex = 0.75
  )
}
