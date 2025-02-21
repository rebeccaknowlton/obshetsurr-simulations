estimate.PTE <- function (df.train, df.test, type, numeric_predictors, categorical_predictors, use.actual.control.S) {
  
  # Get linear estimates
  if (type == "linear" | type == "all") {
    # Fit models without surrogate
    m1 <- lm(as.formula(paste("Y ~", paste(c(numeric_predictors, categorical_predictors), collapse = " + "))),
             data = df.train[df.train$G == 1, ])
    m0 <- lm(as.formula(paste("Y ~", paste(c(numeric_predictors, categorical_predictors), collapse = " + "))),
             data = df.train[df.train$G == 0, ])
    
    # Calculate delta
    delta <- predict(m1, newdata = df.test) - predict(m0, newdata = df.test)
    
    # Fit models including surrogate
    m1.s <- lm(as.formula(paste("Y ~", paste(c(numeric_predictors, categorical_predictors, "S"), collapse = " + "))),
               data = df.train[df.train$G == 1, ])
    m0.s <- lm(as.formula(paste("Y ~", paste(c(numeric_predictors, categorical_predictors, "S"), collapse = " + "))),
               data = df.train[df.train$G == 0, ])
    s0 <- lm(as.formula(paste("S ~", paste(c(numeric_predictors, categorical_predictors), collapse = " + "))),
             data = df.train[df.train$G == 0, ])
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = df.test)
    new.data <- df.test 
    new.data$S <- s.preds
    if (use.actual.control.S) {
      # if user specifies, use the observed s0 values where available
      new.data$S[new.data$G==0] <- df.test$S[df.test$G==0]
    }
    
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
    m1 <- gam(as.formula(paste("Y ~", paste(c(paste0("s(", numeric_predictors, ")"), categorical_predictors), collapse = " + "))), 
              data = df.train[df.train$G == 1, ])
    m0 <- gam(as.formula(paste("Y ~", paste(c(paste0("s(", numeric_predictors, ")"), categorical_predictors), collapse = " + "))),
              data = df.train[df.train$G == 0, ])
    
    # Calculate delta
    delta <- predict(m1, newdata = df.test) - predict(m0, newdata = df.test)
    
    # Fit models including surrogate
    m1.s <- gam(as.formula(paste("Y ~", paste(c(paste0("s(", numeric_predictors, ")"), "s(S)", categorical_predictors), collapse = " + "))), 
                data = df.train[df.train$G == 1, ])
    m0.s <- gam(as.formula(paste("Y ~", paste(c(paste0("s(", numeric_predictors, ")"), "s(S)", categorical_predictors), collapse = " + "))), 
                data = df.train[df.train$G == 0, ])
    s0 <- gam(as.formula(paste("S ~", paste(c(paste0("s(", numeric_predictors, ")"), categorical_predictors), collapse = " + "))), 
              data = df.train[df.train$G == 0, ])
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = df.test)
    new.data <- df.test 
    new.data$S <- s.preds
    if (use.actual.control.S) {
      # if user specifies, use the observed s0 values where available
      new.data$S[new.data$G==0] <- df.test$S[df.test$G==0]
    }
    
    # Calculate delta.s
    delta.s <- predict(m1.s, newdata = new.data) - predict(m0.s, newdata = new.data)
    
    if (type == "all") {
      delta.gam <- delta
      delta.s.gam <- delta.s
    }
  }
  
  # get tree estimates
  if (type == "trees" | type == "all") {
    
    if (length(categorical_predictors) > 0) {
      # One-hot encode categorical variables
      train_cats <- model.matrix(~ . - 1, data = as.data.frame(df.train[, categorical_predictors]))
      test_cats <- model.matrix(~ . - 1, data = as.data.frame(df.test[,categorical_predictors]))
      
      # Combine numeric and encoded categorical predictors
      X_train <- cbind(as.matrix(df.train[, numeric_predictors]), train_cats)
      X_test <- cbind(as.matrix(df.test[, numeric_predictors]), test_cats)
    } else {
      X_train <- as.matrix(df.train[, numeric_predictors])
      X_test <- as.matrix(df.test[, numeric_predictors])
    }
    
    # Fit models without surrogate
    m1 <- regression_forest(
      X = X_train[df.train$G == 1, ],
      Y = df.train[df.train$G == 1, "Y"]
    )
    m0 <- regression_forest(
      X = X_train[df.train$G == 0, ],
      Y = df.train[df.train$G == 0, "Y"]
    )
    
    # Calculate delta
    delta <- predict(m1, newdata = X_test)$predictions -
      predict(m0, newdata = X_test)$predictions
    
    # Fit models including surrogate
    X_train_s <- cbind(X_train, S = df.train$S)
    X_test_s <- cbind(X_test, S = df.test$S)
    
    m1.s <- regression_forest(
      X = X_train_s[df.train$G == 1, ],
      Y = df.train[df.train$G == 1, "Y"]
    )
    
    m0.s <- regression_forest(
      X = X_train_s[df.train$G == 0, ],
      Y = df.train[df.train$G == 0, "Y"]
    )
    
    # Fit surrogate model
    s0 <- regression_forest(
      X = X_train[df.train$G == 0, ],
      Y = df.train[df.train$G == 0, "S"]
    )
    
    # Predict surrogate S in control group using test data
    s.preds <- predict(s0, newdata = X_test)$predictions
    new.data <- X_test_s
    new.data[, "S"] <- s.preds
    if (use.actual.control.S) {
      # if user specifies, use the observed s0 values where available
      new.data[df.test$G==0,"S"] <- df.test$S[df.test$G==0]
    }
    
    # Calculate delta.s
    delta.s <- predict(m1.s, newdata = new.data)$predictions -
      predict(m0.s, newdata = new.data)$predictions
    
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
