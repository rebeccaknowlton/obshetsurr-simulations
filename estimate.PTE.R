estimate.PTE <- function (df.train, df.test, type, predictor_columns) {
  # Get linear estimates
  if (type == "linear" | type == "all") {
    # Fit models without surrogate
    m1 <- lm(as.formula(paste("Y ~", paste(c(predictor_columns), collapse = " + "))),
             data = df.train[df.train$G == 1, ])
    m0 <- lm(as.formula(paste("Y ~", paste(c(predictor_columns), collapse = " + "))),
             data = df.train[df.train$G == 0, ])
    
    # Calculate delta
    delta <- predict(m1, newdata = df.test) - predict(m0, newdata = df.test)
    
    # Fit models including surrogate
    m1.s <- lm(as.formula(paste("Y ~", paste(c(predictor_columns, "S"), collapse = " + "))),
               data = df.train[df.train$G == 1, ])
    m0.s <- lm(as.formula(paste("Y ~", paste(c(predictor_columns, "S"), collapse = " + "))),
               data = df.train[df.train$G == 0, ])
    s0 <- lm(as.formula(paste("S ~", paste(predictor_columns, collapse = " + "))),
             data = df.train[df.train$G == 0, ])
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = df.test)
    new.data <- df.test 
    new.data$S <- s.preds
    #new.data <- cbind(df.test[,colnames(df.test!="S")], s.preds)
    #names(new.data)[dim(new.data)[2]] = "S"
    
    # Calculate delta.s
    delta.s <- predict(m1.s, newdata = new.data) - predict(m0.s, newdata = new.data)
    
    if (type == "all") {
      delta.linear <- delta
      delta.s.linear <- delta.s
    }
  }
  
  # get gam estimates
  if (type == "gam" | type == "all") {
    # Fit models without surrogate
    m1 <- gam(as.formula(paste("Y ~", paste("s(", predictor_columns, ")", collapse = " + "))),
              data = df.train[df.train$G == 1, ])
    m0 <- gam(as.formula(paste("Y ~", paste("s(", predictor_columns, ")", collapse = " + "))),
              data = df.train[df.train$G == 0, ])
    
    # Calculate delta
    delta <- predict(m1, newdata = df.test) - predict(m0, newdata = df.test)
    
    # Fit models including surrogate
    m1.s <- gam(as.formula(paste("Y ~", paste(c(paste("s(", predictor_columns, ")"), "s(S)"), collapse = " + "))),
                data = df.train[df.train$G == 1, ])
    m0.s <- gam(as.formula(paste("Y ~", paste(c(paste("s(", predictor_columns, ")"), "s(S)"), collapse = " + "))),
                data = df.train[df.train$G == 0, ])
    s0 <- gam(as.formula(paste("S ~", paste("s(", predictor_columns, ")", collapse = " + "))),
              data = df.train[df.train$G == 0, ])
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = df.test)
    new.data <- df.test 
    new.data$S <- s.preds
    
    # Calculate delta.s
    delta.s <- predict(m1.s, newdata = new.data) - predict(m0.s, newdata = new.data)
    
    if (type == "all") {
      delta.gam <- delta
      delta.s.gam <- delta.s
    }
  }
  
  # get tree estimates
  if (type == "trees" | type == "all") {
    # Fit models without surrogate
    m1 <- regression_forest(
      X = as.matrix(df.train[df.train$G == 1, predictor_columns]),
      Y = df.train[df.train$G == 1, "Y"]
    )
    m0 <- regression_forest(
      X = as.matrix(df.train[df.train$G == 0, predictor_columns]),
      Y = df.train[df.train$G == 0, "Y"]
    )
    
    # Calculate delta
    delta <- predict(m1, newdata = as.matrix(df.test[, predictor_columns]))$predictions -
      predict(m0, newdata = as.matrix(df.test[, predictor_columns]))$predictions
    
    # Fit models including surrogate
    m1.s <- regression_forest(
      X = as.matrix(df.train[df.train$G == 1, c(predictor_columns, "S")]),
      Y = df.train[df.train$G == 1, "Y"]
    )
    
    m0.s <- regression_forest(
      X = as.matrix(df.train[df.train$G == 0, c(predictor_columns, "S")]),
      Y = df.train[df.train$G == 0, "Y"]
    )
    
    # Fit surrogate model
    s0 <- regression_forest(
      X = as.matrix(df.train[df.train$G == 0, predictor_columns]),
      Y = df.train[df.train$G == 0, "S"]
    )
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = as.matrix(df.test[, predictor_columns]))$predictions
    new.data <- df.test 
    new.data$S <- s.preds
    
    # Calculate delta.s
    delta.s <- predict(m1.s, newdata = as.matrix(new.data[, c(predictor_columns, "S")]))$predictions -
      predict(m0.s, newdata = as.matrix(new.data[, c(predictor_columns, "S")]))$predictions
    
    if (type == "all") {
      delta.trees <- delta
      delta.s.trees <- delta.s
    }
  }
  
  # Add columns to test data and return
  if (type %in% c("linear", "gam", "trees")) {
    df.test$delta <- delta
    df.test$delta.s <- delta.s
    df.test$R.s <- 1 - delta.s / delta
  } else {
    # return all estimates
    df.test$delta.linear <- delta.linear
    df.test$delta.s.linear <- delta.s.linear
    df.test$R.s.linear <- 1 - delta.s.linear / delta.linear
    df.test$delta.gam <- delta.gam
    df.test$delta.s.gam <- delta.s.gam
    df.test$R.s.gam <- 1 - delta.s.gam / delta.gam
    df.test$delta.trees <- delta.trees
    df.test$delta.s.trees <- delta.s.trees 
    df.test$R.s.trees <- 1 - delta.s.trees / delta.trees
  }
  return(df.test) 
}
