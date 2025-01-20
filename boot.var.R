boot.var <- function(df.train, df.test, type, predictor_columns, threshold) {
  num.boot <- 200
  
  if (type %in% c("linear", "gam", "trees")) {
    boot.delta <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
  } else {
    # dataframes for all types
    boot.delta.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.linear <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
    boot.delta.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.gam <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
    boot.delta.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.delta.s.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))
    boot.R.s.trees <- data.frame(matrix(nrow = nrow(df.test), ncol = num.boot))  
  }
  
  for (j in 1:num.boot) {
    boot.data.train <- df.train[sample(1:nrow(df.train), nrow(df.train), replace = TRUE),]
    
    boot.results <- estimate.PTE(df.train = boot.data.train, df.test = df.test, type = type, predictor_columns)
    if (type %in% c("linear", "gam", "trees")) {
      boot.delta[,j] <- boot.results$delta
      boot.delta.s[,j] <- boot.results$delta.s
      boot.R.s[,j] <- boot.results$R.s
    } else {
      boot.delta.linear[,j] <- boot.results$delta.linear
      boot.delta.s.linear[,j] <- boot.results$delta.s.linear
      boot.R.s.linear[,j] <- boot.results$R.s.linear
      boot.delta.gam[,j] <- boot.results$delta.gam
      boot.delta.s.gam[,j] <- boot.results$delta.s.gam
      boot.R.s.gam[,j] <- boot.results$R.s.gam
      boot.delta.trees[,j] <- boot.results$delta.trees
      boot.delta.s.trees[,j] <- boot.results$delta.s.trees
      boot.R.s.trees[,j] <- boot.results$R.s.trees
    }
  }
  
  if (type %in% c("linear", "gam", "trees")) {
    df.test$delta.var <- rowVars(as.matrix(boot.delta))
    df.test$delta.s.var <- rowVars(as.matrix(boot.delta.s))
    df.test$R.s.var <- rowVars(as.matrix(boot.R.s))
  } else {
    df.test$delta.var.linear <- rowVars(as.matrix(boot.delta.linear))
    df.test$delta.s.var.linear <- rowVars(as.matrix(boot.delta.s.linear))
    df.test$R.s.var.linear <- rowVars(as.matrix(boot.R.s.linear))
    df.test$delta.var.gam <- rowVars(as.matrix(boot.delta.gam))
    df.test$delta.s.var.gam <- rowVars(as.matrix(boot.delta.s.gam))
    df.test$R.s.var.gam <- rowVars(as.matrix(boot.R.s.gam))
    df.test$delta.var.trees <- rowVars(as.matrix(boot.delta.trees))
    df.test$delta.s.var.trees <- rowVars(as.matrix(boot.delta.s.trees))
    df.test$R.s.var.trees <- rowVars(as.matrix(boot.R.s.trees))
  }
  
  if (!is.null(threshold)) {
    if (type %in% c("linear", "gam", "trees")) {
      df.test$p.val <- p.adjust(rowMeans(boot.R.s > threshold), method = "BH")
    } else {
      df.test$p.val.linear <- p.adjust(rowMeans(boot.R.s.linear > threshold), method = "BH")
      df.test$p.val.gam <- p.adjust(rowMeans(boot.R.s.gam > threshold), method = "BH")
      df.test$p.val.trees <- p.adjust(rowMeans(boot.R.s.trees > threshold), method = "BH")
    }
  }
  
  return(df.test)
}
