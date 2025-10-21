############################################################
# Section 2.2: Model Fitting Utilities
# Description:
#   Requires "train_data" data.frame
#   Fits a selection of classification models (GLM, LDA, SVM, GBM, rpart,
#   GAM, NNET) to the provided training data, based on boolean flags.
#   The chosen model object and its identifier are returned.
#
# Notes:
#   * Assumes the following objects exist in the environment:
#       - var_names: character vector of feature names (length = num_vars)
#       - num_vars: integer, number of feature variables
#       - lda, svm_u, gbm_used, rpart_used, gam_used, nnet_used: logical flags
############################################################


#############################
# 1.1) SVM Probability Wrapper
#############################
#' Predict class-1 probabilities using an e1071 SVM model
#'
#' @param model_use Model object (e1071::svm) trained with probability=TRUE
#' @param data Newdata for prediction (data.frame)
#'
#' @return Numeric vector of class-1 probabilities
svm_pred <- function(model_use = model_use, data){
  prob_heading <- colnames(attr(predict(model_use, newdata = data, probability = TRUE), "probabilities"))
  prob_nr <- which(prob_heading == "1")[1]
  preds <- attr(predict(model_use, newdata = data, probability = TRUE), "probabilities")[, prob_nr]
  return(preds)
}


########################################
# 1.2) Generic Prediction Function Wrapper
########################################
#' Standardize prediction calls across model families
#'
#' Applies rounding and supports different model families.
#'
#' @param model_use Model object
#' @param model_select Character flag identifying model family
#'   (e.g., "NNET", "RForrest", "LDA", "SVM").
#' @param rounding Integer number of decimal places (floor-based rounding)
#' @param new_data Data frame for prediction
#'
#' @return Numeric vector of predicted probabilities/scores
pred_function <-  function(model_use,
                           model_select,
                           rounding = 2,
                           new_data) {
  
  m_options <- c("NNET", "RForrest", "LDA", "SVM")
  
  if(!(model_select %in% c(m_options)))
    pred_value <- floor(predict(model_use, newdata = new_data, type = "response") * (10^rounding)) / (10^rounding)
  
  if(model_select == "NNET")
    pred_value <- floor(predict(model_use, newdata = new_data) * (10^rounding)) / (10^rounding)
  
  if(model_select == "RForrest")
    pred_value <- floor(predict(model_use, newdata = new_data)[, 2] * (10^rounding)) / (10^rounding)
  
  if(model_select == "LDA")
    pred_value <- floor(predict(model_use, newdata = new_data)$posterior[, 2] * (10^rounding)) / (10^rounding)
  
  if(model_select == "SVM")
    pred_value <- floor(svm_pred(model_use, new_data) * (10^rounding)) / (10^rounding)
  
  return(pred_value)
}



#############################
# 2.1) Modeling Hints (Optional)
#############################
# You can specify the model directly and skip the generic parts by setting:
#   model_select <- "<MODEL>"  # e.g., "GLM", "LDA", "SVM", etc.
#   model_use    <- model_<model>


#################################
# 2.2) Model Fitting Wrapper
#################################
#' Fit a classification model according to selected flags
#'
#' Builds a formula class ~ x1 + x2 + ... and fits the chosen model(s).
#' The last enabled model flag determines the returned model.
#'
#' @param train_data Data frame containing features in columns var_names and
#'   a numeric binary 'class' column.
#'
#' @return list(model_use = <model object>, model_select = <character label>)
model_fitting <- function(train_data = train_data) {

  # Build formula: class ~ x1 + x2 + ...
  formula_data <- paste("class ~ ", var_names[1])
  for (i in 2:num_vars) {
    formula_data <- paste(formula_data, " + ", var_names[i])
  }
  formula_data <- as.formula(formula_data)

  # Base model: GLM (logistic)
  if (glm_use == T) {  
    model_select <- "GLM"
    model_glm <- glm(formula_data, data = train_data, family = binomial(link = "logit"))
    model_use <- model_glm
  }

  # LDA (if enabled)
  if (lda == T) {
    model_lda <- lda(formula_data, data = train_data)
    model_select <- "LDA"
    # model_qda <- qda(formula_data, data = train_data)
  }

  # Other models when LDA flag is FALSE
  if (lda == F) {

    # SVM (if enabled)
    # Example alt options (commented): kernel='linear', scale=FALSE
    if (svm_u == T)  {
      model_svm <- svm(formula_data, data = train_data, probability = TRUE, type = "C-classification")
      model_select <- "SVM"
    }

    # GBM (if enabled)
    if (gbm_used == T) {
      model_gbm <- gbm(
        formula_data,
        data = train_data,
        distribution = "bernoulli",
        shrinkage = 0.15,
        n.minobsinnode = 10,
        n.trees = 500,
        interaction.depth = 3
      )
      model_select <- "GBM"
    }

    # rpart (if enabled)
    if (rpart_used == T) {
      train_data_rpart <- train_data
      train_data_rpart[, num_vars + 1] <- as.factor(train_data_rpart[, num_vars + 1])
      model_rpart <- rpart(formula_data, data = train_data, method = "class", model = T, minsplit = 5)
      model_select <- "RForrest"
    }

    # GAM (if enabled)
    if (gam_used == T) {
      formula_data_gam <- paste("class ~ s(", var_names[1], ")")
      for (i in 2:num_vars) {
        formula_data_gam <- paste(formula_data_gam, " + s(", var_names[i], ")")
      }
      formula_data_gam <- as.formula(formula_data_gam)
      model_gam <- gam(formula_data_gam, data = train_data, family = binomial(link = "logit"))
      model_select <- "GAM"
    }

    # NNET (if enabled)
    if (nnet_used == T) {
      library(nnet)
      model_nn <- nnet(
        formula_data,
        data = train_data,
        size = c(20),
        decay = 0.001,
        maxit = 1000,
        linout = FALSE,  # Binary classification (no linear output)
        trace = FALSE
      )
      model_select <- "NNET"
    }
  }

  # Select final model object according to flags
  if (lda == T)              model_use <- model_lda   # or model_qda
  if (gam_used == T)         model_use <- model_gam
  if (gbm_used == T)         model_use <- model_gbm
  if (rpart_used == T)       model_use <- model_rpart
  if (svm_u == T)            model_use <- model_svm
  if (nnet_used == T)        model_use <- model_nn

  return(list(
    model_use   = model_use,
    model_select = model_select
  ))
}

#########################
####### END #############
#########################