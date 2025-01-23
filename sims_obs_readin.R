# run all the functions in the master file
# set the working directory to where the txt files are
# make sure the parameters below are up-to-date with the simulations

setwd("C:/Users/rkkno/Documents/University of Texas at Austin/Observational/output files")
library(quantreg)
library(ggplot2)
setting <- 1
k <- 0.5
n.test <- 200  



outputfile=c()
for(u in 1:20) {
  outputfile = rbind(outputfile,read.table(paste("obs_outputfile", setting, "_011725_",u,".txt", sep=""), header = T))
}

n.iter <- nrow(outputfile) / n.test

# Generate X1.grid and calculate truth
set.seed(1)
grid.size <- 20
X1.grid.large <- runif(10000, min = get.parameters(setting)$X.distributions$X1$min, 
                       max = get.parameters(setting)$X.distributions$X1$max)
X1.grid <- seq(from = quantile(X1.grid.large, 0.1), to = quantile(X1.grid.large, 0.9), length = grid.size)
truth <- get.truth(setting, X1.grid)

# Matrices to store interpolated results
delta.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.var.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.var.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.var.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.lower.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.lower.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.lower.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.upper.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.upper.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.upper.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)

delta.s.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.var.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.var.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.var.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.lower.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.lower.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.lower.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.upper.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.upper.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
delta.s.upper.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)

R.s.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.var.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.var.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.var.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.lower.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.lower.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.lower.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.upper.grid.lin <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.upper.grid.gam <- matrix(NA, nrow = n.iter, ncol = grid.size)
R.s.upper.grid.tree <- matrix(NA, nrow = n.iter, ncol = grid.size)

# Interpolate each iteration on X1.grid
for (i in 1:n.iter) {
  idx <- ((i - 1) * n.test + 1):(i * n.test)  # Get rows for this iteration
  X1.test <- outputfile$X1[idx]
  
  # Interpolate estimates
  delta.grid.lin[i, ] <- approx(X1.test, outputfile$delta.linear[idx], xout = X1.grid)$y
  delta.grid.gam[i, ] <- approx(X1.test, outputfile$delta.gam[idx], xout = X1.grid)$y
  delta.grid.tree[i, ] <- approx(X1.test, outputfile$delta.trees[idx], xout = X1.grid)$y
  delta.var.grid.lin[i, ] <- approx(X1.test, outputfile$delta.var.linear[idx], xout = X1.grid)$y
  delta.var.grid.gam[i, ] <- approx(X1.test, outputfile$delta.var.gam[idx], xout = X1.grid)$y
  delta.var.grid.tree[i, ] <- approx(X1.test, outputfile$delta.var.trees[idx], xout = X1.grid)$y
  delta.lower.grid.lin[i, ] <- approx(X1.test, outputfile$delta.lower.linear[idx], xout = X1.grid)$y
  delta.lower.grid.gam[i, ] <- approx(X1.test, outputfile$delta.lower.gam[idx], xout = X1.grid)$y
  delta.lower.grid.tree[i, ] <- approx(X1.test, outputfile$delta.lower.trees[idx], xout = X1.grid)$y
  delta.upper.grid.lin[i, ] <- approx(X1.test, outputfile$delta.upper.linear[idx], xout = X1.grid)$y
  delta.upper.grid.gam[i, ] <- approx(X1.test, outputfile$delta.upper.gam[idx], xout = X1.grid)$y
  delta.upper.grid.tree[i, ] <- approx(X1.test, outputfile$delta.upper.trees[idx], xout = X1.grid)$y
  
  delta.s.grid.lin[i, ] <- approx(X1.test, outputfile$delta.s.linear[idx], xout = X1.grid)$y
  delta.s.grid.gam[i, ] <- approx(X1.test, outputfile$delta.s.gam[idx], xout = X1.grid)$y
  delta.s.grid.tree[i, ] <- approx(X1.test, outputfile$delta.s.trees[idx], xout = X1.grid)$y
  delta.s.var.grid.lin[i, ] <- approx(X1.test, outputfile$delta.s.var.linear[idx], xout = X1.grid)$y
  delta.s.var.grid.gam[i, ] <- approx(X1.test, outputfile$delta.s.var.gam[idx], xout = X1.grid)$y
  delta.s.var.grid.tree[i, ] <- approx(X1.test, outputfile$delta.s.var.trees[idx], xout = X1.grid)$y
  delta.s.lower.grid.lin[i, ] <- approx(X1.test, outputfile$delta.s.lower.linear[idx], xout = X1.grid)$y
  delta.s.lower.grid.gam[i, ] <- approx(X1.test, outputfile$delta.s.lower.gam[idx], xout = X1.grid)$y
  delta.s.lower.grid.tree[i, ] <- approx(X1.test, outputfile$delta.s.lower.trees[idx], xout = X1.grid)$y
  delta.s.upper.grid.lin[i, ] <- approx(X1.test, outputfile$delta.s.upper.linear[idx], xout = X1.grid)$y
  delta.s.upper.grid.gam[i, ] <- approx(X1.test, outputfile$delta.s.upper.gam[idx], xout = X1.grid)$y
  delta.s.upper.grid.tree[i, ] <- approx(X1.test, outputfile$delta.s.upper.trees[idx], xout = X1.grid)$y
  
  R.s.grid.lin[i, ] <- approx(X1.test, outputfile$R.s.linear[idx], xout = X1.grid)$y
  R.s.grid.gam[i, ] <- approx(X1.test, outputfile$R.s.gam[idx], xout = X1.grid)$y
  R.s.grid.tree[i, ] <- approx(X1.test, outputfile$R.s.trees[idx], xout = X1.grid)$y
  R.s.var.grid.lin[i, ] <- approx(X1.test, outputfile$R.s.var.linear[idx], xout = X1.grid)$y
  R.s.var.grid.gam[i, ] <- approx(X1.test, outputfile$R.s.var.gam[idx], xout = X1.grid)$y
  R.s.var.grid.tree[i, ] <- approx(X1.test, outputfile$R.s.var.trees[idx], xout = X1.grid)$y
  R.s.lower.grid.lin[i, ] <- approx(X1.test, outputfile$R.s.lower.linear[idx], xout = X1.grid)$y
  R.s.lower.grid.gam[i, ] <- approx(X1.test, outputfile$R.s.lower.gam[idx], xout = X1.grid)$y
  R.s.lower.grid.tree[i, ] <- approx(X1.test, outputfile$R.s.lower.trees[idx], xout = X1.grid)$y
  R.s.upper.grid.lin[i, ] <- approx(X1.test, outputfile$R.s.upper.linear[idx], xout = X1.grid)$y
  R.s.upper.grid.gam[i, ] <- approx(X1.test, outputfile$R.s.upper.gam[idx], xout = X1.grid)$y
  R.s.upper.grid.tree[i, ] <- approx(X1.test, outputfile$R.s.upper.trees[idx], xout = X1.grid)$y
}

# Function to calculate statistics for table
calculate.summaries <- function(estimate.grid, variance.grid, lower.grid, upper.grid, truth.vec) {
  
  truth.mat <- matrix(rep(truth.vec, each = n.iter), nrow = n.iter)
  ESE <- apply(estimate.grid, 2, sd)
  Bias <- colMeans(estimate.grid - truth.mat)
  ASE <- colMeans(sqrt(variance.grid))
  MSE <- colMeans((estimate.grid - truth.mat)^2) 
  Coverage <- colMeans((lower.grid <= truth.mat) & (upper.grid >= truth.mat))
  CI.lower <- colMeans(lower.grid)
  CI.upper <- colMeans(upper.grid)
  #Coverage <- colMeans((estimate.grid - 1.96 * sqrt(variance.grid) <= truth.mat) & 
  #                       (estimate.grid + 1.96 * sqrt(variance.grid) >= truth.mat))
  
  # Combine metrics
  result <- rbind(Bias, ESE, ASE, MSE, Coverage, CI.lower, CI.upper)
  result <- cbind(result, rowMeans(result))  # Add average over X1.grid
  colnames(result) <- c(as.character(round(X1.grid, 2)), "Avg")
  return(result)
}

# Create results tables
results.delta.linear <- calculate.summaries(delta.grid.lin, delta.var.grid.lin, delta.lower.grid.lin, delta.upper.grid.lin, truth$delta)
results.delta.gam <- calculate.summaries(delta.grid.gam, delta.var.grid.gam, delta.lower.grid.gam, delta.upper.grid.gam, truth$delta)
results.delta.tree <- calculate.summaries(delta.grid.tree, delta.var.grid.tree, delta.lower.grid.tree, delta.upper.grid.tree, truth$delta)
results.delta.s.linear <- calculate.summaries(delta.s.grid.lin, delta.s.var.grid.lin, delta.s.lower.grid.lin, delta.s.grid.lin, truth$delta.s)
results.delta.s.gam <- calculate.summaries(delta.s.grid.gam, delta.s.var.grid.gam, delta.s.lower.grid.gam, delta.s.upper.grid.gam, truth$delta.s)
results.delta.s.tree <- calculate.summaries(delta.s.grid.tree, delta.s.var.grid.tree, delta.s.lower.grid.tree, delta.s.upper.grid.tree, truth$delta.s)
results.R.s.linear <- calculate.summaries(R.s.grid.lin, R.s.var.grid.lin, R.s.lower.grid.lin, R.s.upper.grid.lin, truth$R.s)
results.R.s.gam <- calculate.summaries(R.s.grid.gam, R.s.var.grid.gam, R.s.lower.grid.gam, R.s.upper.grid.gam, truth$R.s)
results.R.s.tree <- calculate.summaries(R.s.grid.tree, R.s.var.grid.tree, R.s.lower.grid.tree, R.s.upper.grid.tree, truth$R.s)


# select indices for smaller table: 2 boundary, 2 closer to middle
indices <- c(2,8,12,19, 21)
results.delta.linear.small <- results.delta.linear[,indices]
results.delta.gam.small <- results.delta.gam[,indices]
results.delta.tree.small <- results.delta.tree[,indices]
results.delta.s.linear.small <- results.delta.s.linear[,indices]
results.delta.s.gam.small <- results.delta.s.gam[,indices]
results.delta.s.tree.small <- results.delta.s.tree[,indices]
results.R.s.linear.small <- results.R.s.linear[,indices]
results.R.s.gam.small <- results.R.s.gam[,indices]
results.R.s.tree.small <- results.R.s.tree[,indices]
  
# make combined table with averages only and drop CI lower/upper
table.all <- data.frame(
  Linear = results.R.s.linear[1:5, "Avg"],
  GAM = results.R.s.gam[1:5, "Avg"],
  Trees = results.R.s.tree[1:5, "Avg"]
)

# Print the summary table
print(round(table.all, 3))
latex.table(format(round(as.matrix(table.all),3), nsmall =3), paste0("simres_",setting,"_R.s"), caption = "", dcolumn = T)

# ID strong surrogate
cutoff <- get.X1.cutoff(setting, k)
outputfile$true.flag <- outputfile$X1 >= cutoff  

# Function to compute sensitivity, specificity, PPV, and NPV
flag.summaries <- function(pvals, true.flags) {
  predicted_flag <- (pvals < 0.05)  # Reject null if p-value < 0.05
  
  TP <- sum(predicted_flag & true.flags)   # True Positives
  FP <- sum(predicted_flag & !true.flags)  # False Positives
  TN <- sum(!predicted_flag & !true.flags) # True Negatives
  FN <- sum(!predicted_flag & true.flags)  # False Negatives
  
  sensitivity <- TP / (TP + FN)  # True Positive Rate
  specificity <- TN / (TN + FP)  # True Negative Rate
  PPV <- TP / (TP + FP)  # Precision
  NPV <- TN / (TN + FN)  # Negative Predictive Value
  
  return(c(Sensitivity = sensitivity, Specificity = specificity, 
           PPV = PPV, NPV = NPV))
}

flag.results <- cbind(flag.summaries(outputfile$p.val.linear, outputfile$true.flag),
                      flag.summaries(outputfile$p.val.gam, outputfile$true.flag),
                      flag.summaries(outputfile$p.val.trees, outputfile$true.flag))
colnames(flag.results) <- c("Linear", "GAM", "Trees")

# Print the summary table
print(round(flag.results, 3))
latex.table(format(round(as.matrix(flag.results),3), nsmall =3), paste0("IDregion_res_",setting,""), caption = "", dcolumn = T)


# Function to prepare data for plotting
prepare_plot_data <- function(results, method_name) {
  df <- as.data.frame(t(results[,-(grid.size+1)]))  # Transpose to make X1.grid the first column and drop avg col
  df$X1.grid <- as.numeric(rownames(df))  # Convert row names to numeric
  df$Method <- method_name  # Add method identifier
  return(df)
}

# Combine only R.s results into a single dataframe for plotting
plot_data_Rs <- rbind(
  prepare_plot_data(results.R.s.linear, "Linear"),
  prepare_plot_data(results.R.s.gam, "GAM"),
  prepare_plot_data(results.R.s.tree, "Trees")
)
plot_data_Rs$Method <- factor(plot_data_Rs$Method, levels = c("Linear", "GAM", "Trees"))
plot_data_Rs$Truth <- truth$R.s
plot_data_Rs$Estimate <- plot_data_Rs$Truth + plot_data_Rs$Bias
plot_data_Rs$CI.lower <- plot_data_Rs$CI.lower
plot_data_Rs$CI.upper <- plot_data_Rs$CI.upper
#plot_data_Rs$CI.lower <- plot_data_Rs$Estimate - 1.96*plot_data_Rs$ASE
#plot_data_Rs$CI.upper <- plot_data_Rs$Estimate + 1.96*plot_data_Rs$ASE

ggplot(plot_data_Rs, aes(x = X1.grid, y = Estimate)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.2) +
  geom_line(aes(y = Truth), color = "black", linetype = "dashed", size = 1) +
  facet_wrap(~Method)  +
  labs(title = "Estimated PTE",
       x = "X1",
       y = "R") +
  theme_minimal()

# eventually when I have results for all settings, want to make this into a combined plot where each row is a different setting