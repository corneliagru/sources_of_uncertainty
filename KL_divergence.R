library(MASS)
set.seed(123)

# create necessary directories to save files and plots
mainDir <- getwd()
subDir <- "data"
ifelse(!dir.exists(file.path(mainDir, subDir)), 
       dir.create(file.path(mainDir, subDir)), FALSE)

subDir <- "plots"
ifelse(!dir.exists(file.path(mainDir, subDir)), 
       dir.create(file.path(mainDir, subDir)), FALSE)

# functions ---------------------------------------------------------------

pen.inv <- function(A, method = "ginv", pen.par = 1) {
  # method describes the way of penalization.
  # method = ginv is generalized inverse
  # method = ridge is shrinkage/ ridge penalty
  if (method == "ginv") {
    return(ginv(A))
  }
  if (method == "ridge") {
    return(solve(A + pen.par * diag(dim(A)[1])))
  }
}



init.beta <- function(beta_initialisation = "constant", p) {
  # initialize betas
  # constant = all set to 1
  # decreasing = starts with 1, subrtacts 1/150 per next parameter
  # normal = samples from N(0, 0.5^2)
  if (beta_initialisation == "constant") {
    beta <- matrix(1, p, 1)
  } else if (beta_initialisation == "decreasing") {
    beta <- matrix(1, p, 1) - c(1:200) / 150
  } else if (beta_initialisation == "normal") {
    beta = rnorm(p, sd = 0.5)
  }
  
  # Last 50 betas are set to zero
  beta[150:200] <- 0
  
  return(beta)
}


simulate.KL <- function(beta_initialisation, p = 200, method, pen.par = 0.0001, 
                        n.train = 100, n.test = 1000, sd.y = 0.1,  sd.x = 1) {
  # beta_initialisation for init.beta
  # p is total number of covariates for init beta
  # n.train number of training observations
  # n.test number of test observations
  # sd.y standard deviation of predictor y
  # sd.x standard deviation of covariates x

  
  #setup variables
  beta <- init.beta(beta_initialisation, p)
  
  KL <- c()
  KL.1 <- c()
  KL.2 <- c()
  bias <- c()
  
  # To reduce computation time, we simulate every fifth number of parameters
  p.model <- seq(2, 200, by = 5)
  
  print(paste0("simulating with ", beta_initialisation, " coefficients and regularisation via ", method))
  for (pp in p.model) {
    # We include successively more covariates, but not randomly.
    print(paste0("current number of parameters: ", pp))
    pred.i <- c()
    m.i <- c()
    var.i <- c()
    KL.i <- c()
    KL.1.i <- c()
    KL.2.i <- c()
    for (i in 1:100)
    {
      #  Simulating Training Data
      X <- matrix(rnorm(n.train * p, sd = sd.x), nrow = n.train)
      y <- X %*% matrix(beta) + rnorm(n.train, sd = sd.y)
      
      # Index that contains the included covariates
      index <- c(1:pp) 
      
      # Reduce the design matrix to the included covariates
      
      X.p <- X[, index] 
      
      # Regularization:
      # Taking the generalized inverse
      beta.hat <- pen.inv((t(X.p) %*% X.p),
                          method = method, pen.par = pen.par * p) %*% t(X.p) %*% y
      
      # Simulate test data
      X.test <- matrix(rnorm(n.test * p, sd = sd.x), nrow = n.test)
      X.test.p <- X.test[, index]
      mu.test <- X.test %*% matrix(beta)
      y.test <- mu.test + rnorm(n.test, sd = sd.y)
      
      # Best Model based on selected covariates of X
      fit.p.test <- X.test.p %*% beta.hat # fitted prediction
      fit.p.sd <- sqrt(var(y.test - fit.p.test)) # Variance
      
      # Kullback Leibler Divergence KL(f(y|x), p(y |x , \hat{\theta}))
      KL.i <- c(KL.i, sum(log(dnorm(y.test, mu.test, sd = sd.y))) -
                  sum(log(dnorm(y.test, fit.p.test, sd = fit.p.sd))))
      
      # Optimal Parameter
      fit.p.test.0 <-
        X.test.p %*% solve(t(X.test.p) %*% X.test.p) %*% t(X.test.p) %*% mu.test
      fit.p.sd.0 <- sqrt(var(y.test - fit.p.test.0)) # Variance
      
      # part one, relating to "bias", 
      # i.e., "distance" from true relationship to optimal parameters
      KL.1.i <-
        c(KL.1.i, sum(log(dnorm(y.test, mu.test, sd = sd.y))) -
            sum(log(dnorm(y.test, fit.p.test.0, sd = fit.p.sd.0))))
      
      # part two, relating to "variance"
      # i.e., "distance" from optimal parameters to current estimate
      KL.2.i <-
        c(KL.2.i, sum(log(
          dnorm(y.test, fit.p.test.0, sd = fit.p.sd.0))) -
            sum(log(dnorm(y.test, fit.p.test, sd = fit.p.sd))))
    }
    
    KL <- c(KL, mean(KL.i))
    KL.1 <- c(KL.1, mean(KL.1.i))
    KL.2 <- c(KL.2, mean(KL.2.i))
  }
  # save simulation results
  result <- data.frame(p.model = p.model, KL = KL, KL.1 = KL.1, KL.2 = KL.2, 
                       beta_initialisation = beta_initialisation, method = method)
  
  objname <- paste0("KL_", beta_initialisation, "_", method)
  
  assign(objname, result)
  save(list = objname, file = paste0("data/", objname, ".Rdata"))
  print("simulation finished")
  return(result)
}


save_plots <- function(result, hide_yaxt = TRUE, start_at_0 = FALSE, add_to_filename = "") {
  
  objname <- paste0("KL_", result[1, "beta_initialisation"], "_", result[1, "method"], add_to_filename)
  
  pdf(paste0("plots/", objname, ".pdf"),
      width = 5,
      height = 7)
  opar <-
    par(
      mfrow = c(3, 1),
      cex.main = 2,
      cex.lab = 2,
      cex.axis = 2
    )
  if (hide_yaxt) {
    par(yaxt = "n")
  } else {
      par(
        mar = c(5.1,5,4.1,2.1),
        las = 1
      )
  }
  

  plot(
    result$p.model,
    result$KL,
    type = "l",
    main = expression(paste(
      E[X],
      "{KL(f(" %.% "|x), p(" %.% "|x, ",
      hat(theta)[pen],
      "))}, ",
      bold(" Eq. (12)")
    )),
    xlab = "",
    ylab = "",
    ylim = if (start_at_0) c(0, max(result$KL)) else c(min(result$KL), max(result$KL))
  ) 
  plot(
    result$p.model,
    result$KL.1,
    type = "l",
    main = "Eq. (13)",
    xlab = "",
    ylab = ""
  )
  plot(
    result$p.model,
    result$KL.2,
    type = "l",
    main = "Eq. (14)",
    xlab = "p",
    ylab = ""
  )
  dev.off()
  
  par(opar)
  
}


# run simulation ----------------------------------------------------------

KL_decreasing_ridge <- simulate.KL(beta_initialisation = "decreasing", method = "ridge")
KL_decreasing_ginv <- simulate.KL(beta_initialisation = "decreasing", method = "ginv")
KL_constant_ridge <- simulate.KL(beta_initialisation = "constant", method = "ridge")
KL_constant_ginv <- simulate.KL(beta_initialisation = "constant", method = "ginv")



# load results ------------------------------------------------------------

load("data/KL_decreasing_ridge.Rdata")
load("data/KL_decreasing_ginv.Rdata")
load("data/KL_constant_ridge.Rdata")
load("data/KL_constant_ginv.Rdata")


# save plots --------------------------------------------------------------

save_plots(KL_decreasing_ridge)
save_plots(KL_decreasing_ginv)
save_plots(KL_constant_ridge)
save_plots(KL_constant_ginv)


# with y values
save_plots(KL_decreasing_ridge, hide_yaxt = FALSE, add_to_filename = "_showy")
save_plots(KL_decreasing_ginv, hide_yaxt = FALSE, add_to_filename = "_showy")
save_plots(KL_constant_ridge, hide_yaxt = FALSE, add_to_filename = "_showy")
save_plots(KL_constant_ginv, hide_yaxt = FALSE, add_to_filename = "_showy")

# with y axis starting at origin
save_plots(KL_decreasing_ridge, hide_yaxt = FALSE, start_at_0 = TRUE, add_to_filename = "_showy_start0")
save_plots(KL_decreasing_ginv, hide_yaxt = FALSE, start_at_0 = TRUE, add_to_filename = "_showy_start0")
save_plots(KL_constant_ridge, hide_yaxt = FALSE, start_at_0 = TRUE, add_to_filename = "_showy_start0")
save_plots(KL_constant_ginv, hide_yaxt = FALSE, start_at_0 = TRUE, add_to_filename = "_showy_start0")

# sanity check whether KL = KL.1 + KL.2 is true
#plot(p.model, KL.1 + KL.2, type = "l" )