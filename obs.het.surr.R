obs.het.surr <- function(df.train, df.test, type, var.want = FALSE, threshold = NULL) {
  # Check that inputs valid
  if (!is.data.frame(df.train) || !is.data.frame(df.test)) {
    stop("Both df.train and df.test must be data frames.")
  }
  if (!(type %in% c("linear", "gam", "trees", "all"))) {
    stop('The "type" argument must be one of "linear", "gam", "trees", or "all".')
  }
  
  required_columns <- c("G", "S", "Y")
  missing_columns_train <- setdiff(required_columns, colnames(df.train))
  missing_columns_test <- setdiff(required_columns, colnames(df.test))
  
  if (length(missing_columns_train) > 0) {
    stop(paste("df.train is missing the following required columns:", paste(missing_columns_train, collapse = ", ")))
  }
  
  if (length(missing_columns_test) > 0) {
    stop(paste("df.test is missing the following required columns:", paste(missing_columns_test, collapse = ", ")))
  }
  
  predictor_columns <- setdiff(colnames(df.train), required_columns)
  if (length(predictor_columns) == 0) {
    stop("No predictor columns (e.g., X1, X2, ...) found in the data frames.")
  }
  
  missing_predictors_test <- setdiff(predictor_columns, colnames(df.test))
  if (length(missing_predictors_test) > 0) {
    stop(paste(
      "df.test is missing the following predictor columns from df.train:",
      paste(missing_predictors_test, collapse = ", ")
    ))
  }
  
  # calculate point estimates
  df.test <- estimate.PTE(df.train = df.train, df.test = df.test, type = type, predictor_columns = predictor_columns)
  # calculate variance and flag region if desired
  if (var.want == TRUE) {
    df.test <- boot.var(df.train, df.test, type, predictor_columns, threshold)
  }
  return(df.test)
}
