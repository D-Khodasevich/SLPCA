#####################
# SuperLearner PCA Clock Functions
# Dennis Khodasevich


# train_all_clocks: super function for training all 3 clocks
### meth: dna methylation data input
### pheno: phenotype of interest (numeric vector or matrix column)
### sl_library: superlearner library list
### save_data: whether to save individual clock outputs
##### If set to TRUE
###### cpg_clock, pca_clock, sl_pca_clock: names for the cpg, pca, and sl pca clock outputs respectively
train_all_clocks <- function(meth, pheno, sl_library, family_type = "gaussian", lasso_plot=FALSE, 
                             nfold=10, alpha_val=0.5, method = "method.NNLS",  
                             nlambda_val=100, seed=1234, trim_length=15, 
                             save_data = TRUE, 
                             cpg_clock = "CpG_coefs", 
                             pca_clock = "PCA_Clock", 
                             sl_pca_clock = "SL_PCA_Clock"){
  
  starttime <- Sys.time()
  
  meth.na <- sum(is.na(meth))
  pheno.na <- sum(is.na(pheno))
  
  if(meth.na > 0 | pheno.na > 0){
    stop("Error: Check for NAs in Data Inputs. Perform Imputation of Removal of Missing Observations as Necessary")
  }
  
  ################# CpG Clock procedure
  message("Training CpG-based Predictor")
  
  X.bm <- as.big.matrix(meth)
  Y <- pheno
  fit.enet <- cv.biglasso(X.bm, Y, penalty = 'enet', alpha = 0.5, family = family_type, 
                          seed = 1234, nfolds = 10)
  
  if(lasso_plot == TRUE) {
    jpeg("LASSO_Summary.jpeg", quality = '90%')
    par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
    plot(fit.enet, type = "all") # LASSO summary plots
    dev.off()
  }
  
  # subset and save final CpG predictors
  coefs <- as.matrix(coef(fit.enet))
  coefs <- as.data.frame(coefs)
  colnames(coefs) <- c("beta")
  coefs <- coefs %>% 
    filter(beta != 0)
  coefs$cpgs <- rownames(coefs)
  
  if(save_data == TRUE){
    save(coefs, file = paste0(cpg_clock, ".RData"))
  }
  rm(X.bm)
  
  ################# PCA Clock procedure
  message("Training PCA-based Predictor")
  PCA = prcomp(meth,scale.=F)
  TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-trim_length)]
  
  #Train PC clock. Can test different models using different alpha and lambda parameters 
  # (see glmnet documentation)
  cv = cv.glmnet(TrainPCData, pheno, nfolds=nfold, alpha=alpha_val, family=family_type)
  
  # Final model only uses a subset of PCs. Compress your model:
  CalcPCAge <- vector(mode = "list",length = 0)
  temp = as.matrix(coef(cv,s = cv$lambda.min))
  CalcPCAge$model = temp[temp!=0,][-1]
  CalcPCAge$intercept = temp[1,1]
  CalcPCAge$center = PCA$center
  CalcPCAge$rotation = PCA$rotation[,names(CalcPCAge$model)]
  
  if(save_data == TRUE){
    save(CalcPCAge, file = paste0(pca_clock, ".RData"))
  }
  
  ################# PCA SL Clock procedure
  message("Training PCA SuperLearner-based Predictor")
  TrainPCData <- as.data.frame(TrainPCData)
  
  set.seed(seed)
  sl = SuperLearner(Y = pheno, X = TrainPCData, family = family_type,
                    SL.library = sl_library)
  
  sl_pca <- list(sl, PCA)
  
  if(save_data == TRUE){
    save(sl_pca, file = paste0(sl_pca_clock, ".RData"))
  }
  
  endtime <- Sys.time()
  timediff <- endtime - starttime
  
  all_clocks <- list(coefs, CalcPCAge, sl_pca, timediff)
  names(all_clocks) <- c("CpG_Clock", "PCA_Clock", "SL_PCA_Clock", "Runtime")
  return(all_clocks)
  
}

# predict_all_clocks: super function for obtaining predictions from all 3 clocks
### meth: dna methylation data input
### cpg_coefs: cpg clock input
### pca_model: pca clock input
### pca_sl_model: pcal sl clock input
### age_transform: whether to inverse transform age predictions (set Horvath for childhood clock)
predict_all_clocks <- function(cpg_coefs, pca_model, pca_sl_model, meth, family_type = "gaussian", 
                               trim_length=15, age_transform){
  
  # input checks
  cpg_list <- names(pca_model[["center"]])
  meth <- meth[, c(colnames(meth) %in% cpg_list)]
  present_cpgs <- colnames(meth)
  missing <- setdiff(cpg_list, present_cpgs)
  if(length(missing) > 0) {
    message("Not all CpGs are present.", length(missing), " missing. Proceeding with analysis.")
    missing_fillin <- matrix(NA, ncol = length(missing), nrow = nrow(meth))
    colnames(missing_fillin) <- missing
    rownames(missing_fillin) <- rownames(meth)
    meth <- cbind(meth, missing_fillin)
  }
  
  meth <- meth[ , cpg_list]
  
  meth.na <- sum(is.na(meth))
  if(meth.na > 0){
    message("Missing Data: performing mean imputation to fill in missing data")
    meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
    meth <- apply(meth,1,meanimpute)
    meth <- t(meth)
  }
  
  ################## CpG Prediction
  if(length(cpg_coefs) > 2){
    stop("Error: Model coefficient input must only contain a beta coefficient column 
         and CpG ID column named 'cpgs'")
  }
  
  coef_list <- cpg_coefs$cpgs
  int_val <- cpg_coefs[1, "beta"]
  coef_list <- tail(coef_list, n=length(coef_list)-1)
  testing <- as.data.frame(meth) 
  testing <- testing %>% 
    dplyr::select(all_of(coef_list))
  testing <- as.data.frame(t(testing))
  testing$cpgs <- rownames(testing)
  testdat <- left_join(cpg_coefs, testing, by="cpgs")
  testdat <- testdat[-1, ]
  
  testb <- testdat[,1]
  testx <- testdat[,3:length(testdat)]
  testval <- testb*testx
  testval <- colSums(testval)
  testval <- as.data.frame(testval)
  
  testval$testval <- int_val + testval$testval
  
  preds <- testval$testval
  
  if(age_transform == "Horvath") {
    preds <- (exp(preds + (log(22))))-1
  }
  
  ################## PCA Prediction
  PCAge <- sweep(as.matrix(meth),2,pca_model$center) %*% 
    pca_model$rotation %*% pca_model$model + pca_model$intercept
  
  if(age_transform == "Horvath") {
    PCAge <- (exp(PCAge + (log(22))))-1
  }
  
  preds <- cbind(preds, PCAge)
  
  ################## PCA SL Prediction  
  SL <- pca_sl_model[[1]]
  PCA <- pca_sl_model[[2]]
  
  TestPCData = predict(PCA,meth)[,1:(dim(PCA$x)[2]-trim_length)]
  TestPCData <- as.data.frame(TestPCData)
  
  slpred = predict(SL, TestPCData, onlySL = TRUE)
  PCAge_SL <- slpred$pred
  
  if(age_transform == "Horvath") {
    PCAge_SL <- (exp(PCAge_SL + (log(22))))-1    
  }

  preds <- cbind(preds, PCAge_SL)
  colnames(preds) <- c("CpG_Clock", "PCA_Clock", "PCA_SL_Clock")
  preds <- as.data.frame(preds)
  
  return(preds)
}


# glmnet parameters
SL.glmnet05 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.05, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet10 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.10, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet15 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.15, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet25 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.25, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet50 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.50, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet75 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.75, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}


SL.glmnet60 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.60, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet40 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.40, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet90 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.90, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet0 <- function (Y, X, newX, family, obsWeights, id, alpha = 0, nfolds = 10, 
                        nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet100 <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
                          nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  #    .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

comp_plot <- function(X, Y, x_lab, y_lab) {
  
  corr = round_fix(cor(X, Y), digits=3)
  maeval = round_fix(mdae(X, Y), digits=3)
  
  tempdat <- as.data.frame(cbind(X, Y))
  colnames(tempdat) <- c("X", "Y")
  
  a <- ggplot(tempdat, aes(x=X, y=Y)) + 
    geom_point(size=1.5, alpha=0.8) + 
    geom_abline(slope=1) + 
    theme_bw() + 
    xlab(x_lab) + ylab(y_lab) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste0("Corr. = ", corr, "\n","MAE = ", maeval))
  
  return(a)
}

auc_plots <- function(cpg, pca, sl) {
  
  y <- cpg$sensitivities
  x <- 1 - (cpg$specificities)
  cpg_auc <- round(pROC::auc(cpg), digits=3)
  c <- as.data.frame(cbind(y, x)) 
  
  y <- pca$sensitivities
  x <- 1 - (pca$specificities)
  pca_auc <- round(pROC::auc(pca), digits=3)
  p <- as.data.frame(cbind(y, x) )
  
  y <- sl$sensitivities
  x <- 1 - (sl$specificities)
  sl_auc <- round(pROC::auc(sl), digits=3)
  s <- as.data.frame(cbind(y, x) )
  
  c$Predictor <- paste("CpG") 
  p$Predictor <- paste("PCA") 
  s$Predictor <- paste("SL PCA") 
  plotdat <- rbind(c,p,s)
  
  pl <- ggplot(plotdat, aes(x=x, y=y, color=Predictor)) + 
    geom_line(linewidth=1.5) + 
    theme_classic() + 
    xlab("False-positive rate") + ylab("True-positive rate") + 
    annotate("text", x=0.75, y=0.15, 
             label = paste("AUC Estimates", "\n", 
                           "CpG: ", cpg_auc, "\n",  
                           "PCA: ", pca_auc, "\n",  
                           "SL PCA: ", sl_auc))
  return(pl)
}

# round_fix, round while keeping trailing 0's
### based on function from finalfit R package
round_fix = function(x, digits){
  sprintf.arg = paste0("%.", digits, "f")
  x.rounded = do.call(sprintf, list(sprintf.arg, x)) 
  return(x.rounded)
}
