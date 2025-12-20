# load necessary libraries and source files 
library(mgcv) # for gam
library(grf) # for trees
library(matrixStats)
source("obs.het.surr.R") # main function
source("estimate.PTE.R") # point estimates for delta, delta.s, and R.s
source("boot.var.R") # bootstrap variance


get.parameters <- function(setting) { # linear models correct
  if (setting == 1) {
    list(
      X.distributions = list(
        X1 = list(dist = "uniform", min = 0, max = 3),
        X2 = list(dist = "gamma", shape = 2, scale = 2),
        X3 = list(dist = "uniform", min = 0, max = 5),
        X4 = list(dist = "gamma", shape = 3, scale = 1),
        X5 = list(dist = "uniform", min = 0, max = 2),
        X6 = list(dist = "gamma", shape = 1, scale = 1)
      ),
      G.coefficients = c(0.2, 0.3, 0.5, 0.2, 0.4, 0.1), # Coeffs for X1 to X6
      S.coefficients = c(1.5, 0.2, 0.2, 0.3, 0.1, 0.4, 0.3), # Coeffs for G, X1 to X6
      Y.coefficients = c(1, 2, 0.2, 0.5, 0.2, 0.1, 0.3, 0.4, 2), # Coeffs for G, S, X1 to X6, G*X1
      S.sd = c(0.4, 1.8), # SD for the control and treated groups
      Y.sd = 1
    )
  } else if (setting == 2) { # linear models broken, so GAM should be better
    list(
      X.distributions = list(
        X1 = list(dist = "uniform", min = 0, max = 3),
        X2 = list(dist = "gamma", shape = 2, scale = 2),
        X3 = list(dist = "uniform", min = 0, max = 5),
        X4 = list(dist = "gamma", shape = 3, scale = 1),
        X5 = list(dist = "uniform", min = 0, max = 2),
        X6 = list(dist = "gamma", shape = 1, scale = 1)
      ),
      G.coefficients = c(0.2, 0.3, 0.5, 0.2, 0.4, 0.1), # Coeffs for X1 to X6
      S.coefficients = c(1, 0.2, 0.2, 0.3, 0.1, 0.4, 0.3), # Coeffs for G, X1 to X6
      Y.coefficients = c(1, 2, 1, 1, 1, 1, 1, 1, 1.5), # Coeffs for G, S, X1 to X6, G*X1 in the Y.function
      Y.function = function(Y.coefficients, G, S, X1, X2, X3, X4, X5, X6) {
        Y.coefficients[1] * G + Y.coefficients[2] * S + Y.coefficients[3] * sin(X1) +
          Y.coefficients[4] * cos(X2) + Y.coefficients[5]* X3^2 + Y.coefficients[6] * X4 +
          Y.coefficients[7] * log(X5 + 1) + Y.coefficients[8] * sqrt(X6) + 
          Y.coefficients[9] * (X1^2) * G
      },
      S.sd = c(0.4, 1.8), # SD for the control and treated groups
      Y.sd = 1
    )
  } else if (setting == 3) { # GAM broken
    list(
      X.distributions = list(
        X1 = list(dist = "uniform", min = 0, max = 3),
        X2 = list(dist = "gamma", shape = 2, scale = 2),
        X3 = list(dist = "uniform", min = 0, max = 5),
        X4 = list(dist = "gamma", shape = 3, scale = 1),
        X5 = list(dist = "uniform", min = 0, max = 2),
        X6 = list(dist = "gamma", shape = 1, scale = 1)
      ),
      G.coefficients = c(0.2, 0.3, 0.5, 0.2, 0.4, 0.1), # Coeffs for X1 to X6
      S.coefficients = c(1, 0.2, 0.2, 0.3, 0.1, 0.4, 0.3), # Coeffs for G, X1 to X6
      Y.coefficients = c(1, 2, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 1.5), # Coeffs for G, S, X1 to X6, X1*X5, X2*X3, X4*X6, G*X1 in the Y.function
      #Y.function = function(Y.coefficients, G, S, X1, X2, X3, X4, X5, X6) {
      #  Y.coefficients[1] * G + Y.coefficients[2] * S +
      #    Y.coefficients[9] * exp(-X1^2) +
      #    Y.coefficients[10] * (X2 * sin(X3)) +
      #    Y.coefficients[11] * (X4^2 * log(X6 + 1)) +
      #    Y.coefficients[12] * (X1^2) * G
      #},
      Y.function = function(Y.coefficients, G, S, X1, X2, X3, X4, X5, X6) {
        Y.coefficients[1] * G + Y.coefficients[2] * S +  
          Y.coefficients[9] * X1 * X5^2 + Y.coefficients[10] * log(X2 / X3) + Y.coefficients[11] * sin(X4 + X6) +
          Y.coefficients[12] * (X1^2) * G
      },
      S.sd = c(0.4, 1.8), # SD for the control and treated groups
      Y.sd = 1
    )
  } else if (setting == 4) { # same as setting 1, but with no heterogeneity
    list(
      X.distributions = list(
        X1 = list(dist = "uniform", min = 0, max = 3),
        X2 = list(dist = "gamma", shape = 2, scale = 2),
        X3 = list(dist = "uniform", min = 0, max = 5),
        X4 = list(dist = "gamma", shape = 3, scale = 1),
        X5 = list(dist = "uniform", min = 0, max = 2),
        X6 = list(dist = "gamma", shape = 1, scale = 1)
      ),
      G.coefficients = c(0.2, 0.3, 0.5, 0.2, 0.4, 0.1),
      S.coefficients = c(1, 0.2, 0.2, 0.3, 0.1, 0.4, 0.3),
      Y.coefficients = c(1, 2, 0.2, 0.5, 0.2, 0.1, 0.3, 0.4, 0),
      S.sd = c(0.4, 1.8), # SD for the control and treated groups
      Y.sd = 1
    )
  } else {
    stop("Unsupported setting.")
  }
}


generate.data <- function(n, setting) {
  params <- get.parameters(setting)
  
  # Generate baseline covariates based on specified distributions
  generate.variable <- function(distribution) {
    if (distribution$dist == "gamma") {
      rgamma(n, shape = distribution$shape, scale = distribution$scale)
    } else if (distribution$dist == "uniform") {
      runif(n, min = distribution$min, max = distribution$max)
    } else {
      stop("Unsupported distribution.")
    }
  }
  
  X1 <- generate.variable(params$X.distributions$X1)
  X2 <- generate.variable(params$X.distributions$X2)
  X3 <- generate.variable(params$X.distributions$X3)
  X4 <- generate.variable(params$X.distributions$X4)
  X5 <- generate.variable(params$X.distributions$X5)
  X6 <- generate.variable(params$X.distributions$X6)
  
  # Treatment assignment (observational, depends on all covariates)
  G.prob <- plogis(params$G.coefficients[1] * X1 + params$G.coefficients[2] * X2 +
                     params$G.coefficients[3] * X3 + params$G.coefficients[4] * X4 +
                     params$G.coefficients[5] * X5 + params$G.coefficients[6] * X6)
  G <- rbinom(n, 1, G.prob * 0.5)
  
  # Surrogate marker S

  S <- params$S.coefficients[1] * G +
    params$S.coefficients[2] * X1 + params$S.coefficients[3] * X2 +
    params$S.coefficients[4] * X3 + params$S.coefficients[5] * X4 +
    params$S.coefficients[6] * X5 + params$S.coefficients[7] * X6 +
    rnorm(n, 0, params$S.sd[1] * (1 - G) + params$S.sd[2] * G)

  # Outcome Y
  if (setting == 1 | setting == 4) {
    Y <- params$Y.coefficients[1] * G + params$Y.coefficients[2] * S +
      params$Y.coefficients[3] * X1 + params$Y.coefficients[4] * X2 +
      params$Y.coefficients[5] * X3 + params$Y.coefficients[6] * X4 +
      params$Y.coefficients[7] * X5 + params$Y.coefficients[8] * X6 +
      params$Y.coefficients[9] * X1 * G + rnorm(n, 0, params$Y.sd)
  } else {
    Y <- params$Y.function(params$Y.coefficients, G, S, X1, X2, X3, X4, X5, X6) + rnorm(n, 0, params$Y.sd)
  }
  
  data.frame(X1, X2, X3, X4, X5, X6, G, S, Y)
}

get.truth <- function(setting, grid) {
  params <- get.parameters(setting)
  
  if (setting == 1 | setting == 4) {
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[9]
    
    # Calculate truth at each grid point
    true.delta <- sapply(grid, function(x) {
      GX1.coef * x + G.coef + S.coef * (params$S.coefficients[1])
    })
    true.delta.s <- sapply(grid, function(x) {
      GX1.coef * x + G.coef
    })
  }
  
  if (setting == 2) {
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[9]
    
    # Calculate truth at each grid point
    true.delta <- sapply(grid, function(x) {
      GX1.coef * (x^2) + G.coef + S.coef * (params$S.coefficients[1])
    })
    true.delta.s <- sapply(grid, function(x) {
      GX1.coef * (x^2) + G.coef
    })
  }
  
  if (setting == 3) { 
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[12]
    
    # Calculate truth at each grid point
    true.delta <- sapply(grid, function(x) {
      GX1.coef * (x^2) + G.coef + S.coef * (params$S.coefficients[1])
    })
    true.delta.s <- sapply(grid, function(x) {
      GX1.coef * (x^2) + G.coef
    })
  }
  
  true.R.s <- 1 - true.delta.s / true.delta
  data.frame(X1 = grid, delta = true.delta, delta.s = true.delta.s, R.s = true.R.s)
}

# function to get X1 cutoff where PTE > threshold. Need to fill out for other settings later, basically the inverse of the get.truth
get.X1.cutoff <- function(setting, threshold) {
  params <- get.parameters(setting)
  
  if (setting == 1) {
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[9]
    
    # this could change if setting changes! check algebra
    cutoff <- ((( S.coef * (params$S.coefficients[1])) / threshold) - G.coef -  (S.coef * (params$S.coefficients[1]))) / GX1.coef
  }
  
  if (setting == 2) {
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[9]
    
    # this could change if setting changes! check algebra
    cutoff <- sqrt(((( S.coef * (params$S.coefficients[1])) / threshold) - G.coef -  (S.coef * (params$S.coefficients[1]))) / GX1.coef)
  }
  
  if (setting == 3) {
    # Extract coefficients
    G.coef <- params$Y.coefficients[1]
    S.coef <- params$Y.coefficients[2]
    GX1.coef <- params$Y.coefficients[12]
    
    # this could change if setting changes! check algebra
    cutoff <- sqrt(((( S.coef * (params$S.coefficients[1])) / threshold) - G.coef -  (S.coef * (params$S.coefficients[1]))) / GX1.coef)
  }
  
  return(cutoff)
}


#THIS IS WHAT MAKES EACH PARALLEL VERSION DIFFERENT
set.seed(parallel.num*100)

outputfile = c()

for (i in 1:num.sim) {
  start_time <- Sys.time()
  
  data.temp <- generate.data(n = n, setting = setting)
  
  # split into test and train
  test.idx <- sample(1:n, size = n.test, replace = FALSE)
  data.test <- data.temp[test.idx,]
  data.train <- data.temp[-test.idx,]
  
  output <- obs.het.surr(df.train = data.train, df.test = data.test, type = "all", var.want = TRUE, threshold = k)
  outputfile <- rbind(outputfile,output)
  
  print(i)
  end_time <- Sys.time()
  print(end_time-start_time)
}

write.table(outputfile, paste("obs_outputfile", setting, "_011725_",parallel.num,".txt", sep=""), quote = FALSE, row.names = FALSE)



