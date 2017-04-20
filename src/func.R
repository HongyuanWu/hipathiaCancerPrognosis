# file:
#   func.R
# description:
#   implements general functions and classifiers for the entire study



# libraries ---------------------------------------------------------------



# # generic libraries needed for general purpose
# # from CRAN
# library(knitr) # knitr::knit
# library(devtools) # devtools::session_info
# library(caret) # caret::createMultiFolds
# library(ROCR) # classification performance evaluation
# library(ggplot2) # pretty plot
# library(survival) # survival obj
# # from Bioconductor
# library(survcomp) # survcomp::concordance.index
# library(grid) # multiplot with ggplot

# # import classifier libraries in those predictors being called
# # from CRAN
# library(glmnet) # L1-GLM
# library(LiblineaR) # L1-reg L2-loss SVM
# library(gbm) # GBM
# library(kernlab) # kendall kernel svm and estimate gamma for Gaussian RBF SVM
# library(pcaPP) # fast computation kendall kernel by pcaPP::cor.fk
# library(pamr) # PAM
# library(MASS) # LDA
# library(e1071) # NB, SVM and KNN tuning
# library(class) # KNN predicting
# library(randomForest) # RF

library(methods)
options(stringsAsFactors = FALSE)



# data processing ---------------------------------------------------------



normalizeData <- function(xtr, xtst=NULL, do.center=TRUE, do.scale=TRUE, ...)
{
  # normalizeData() normalizes col-wisely to mean 0 (do.center=T) and sd 1 (do.scale=T)
  # NOTE: xtst is normalized wrt xtr
  
  if(is.vector(xtr)){
    warning("vector features are coerced to a matrix!")
    xtr <- matrix(xtr, ncol=1)
    if(!is.null(xtst)){
      xtst <- matrix(xtst, ncol=1)
    }
  }
  
  # only normalize numeric cols
  numid <- sapply(as.data.frame(xtr), is.numeric)
  t <- scale(xtr[,numid], center=do.center, scale=do.scale)
  xtr[,numid] <- t
  
  if(is.null(xtst)){
    return(xtr)
  } else{
    if(do.center){
      xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:center"), FUN="-")
    }
    if(do.scale){
      xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:scale"), FUN="/")
    }
    return(list(xtr=xtr, xtst=xtst))
  }
}



removeConst <- function(xtr, xtst=NULL, tol2const=1e-6, ...)
{
  numid <- which(sapply(as.data.frame(xtr), is.numeric))
  constid <- apply(xtr[,numid], 2, function(u){
    (abs(mean(u, na.rm = TRUE)) < tol2const && diff(range(u, na.rm = TRUE)) < tol2const) || (diff(range(u, na.rm = TRUE)) / abs(mean(u, na.rm = TRUE)) < tol2const)
  })
  xtr <- xtr[ , numid[!constid], drop = F]
  if(is.null(xtst)){
    return(xtr)
  } else{
    xtst <- xtst[ , numid[!constid], drop = F]
    return(list(xtr=xtr, xtst=xtst))
  }
}



computeKernelMatrix <- function(x, kernel, ...)
{
  # @param x n*p data matrix
  # @param kenrel must be a function!
  # @return n*n kernel matrix
  
  stopifnot(is.function(kernel))
  kmat <- kernlab::kernelMatrix(kernel, as.matrix(x))
  dimnames(kmat) <- list(rownames(x),rownames(x))
  return(kmat)
}



centerScaleKernelMatrix <- function(kmat, train.idx = 1:nrow(kmat), ...)
{
  # @param kmat SYMMETRIC kernel matrix
  # @param train.idx Indices for training set in rows/cols and the others are assumed test set
  # @return centered and scaled kernel matrix with respect to training set in feature space
  
  stopifnot(isSymmetric(kmat))
  n <- nrow(kmat)
  ntr <- length(train.idx)
  
  #centering
  ed <- matrix(0, nrow = n, ncol = n, dimnames = dimnames(kmat))
  ed[train.idx, ] <- 1
  kmcs <- kmat - t(ed)%*%kmat/ntr - kmat%*%ed/ntr + t(ed)%*%kmat%*%ed/(ntr^2)
  #scaling
  dsr <- sqrt(diag(kmcs))
  kmcs <- sweep(
    sweep(kmcs, 1, dsr, FUN='/'),
    2, dsr, FUN='/')
  dimnames(kmcs) <- dimnames(kmat)
  kmcs <- kernlab::as.kernelMatrix(kmcs)
  
  return(kmcs)
}



# validation --------------------------------------------------------------



evaluatePred <- function(pred, ytst, ysurv = NULL, pos.label = tail(names(table(ytst)),1), 
                         savepath = NULL, beta = 1, ...)
{
  # pred results from calling predictorXXX(), ess. a list containing $class and $prob of which class must be factor or character
  # ytst is true (binary or multi) class indicators
  # ysurv is of class Surv containing true surv.time and surv.event, defaulted unavailable
  # pos.label defines positive label(s) in case of computing tpr, fpr, etc, defaulted the last level in ytst
  # NOTE large pred$prob should be associated with short survival time!!
  # NOTE F(beta=1) is used by default while F(beta=0.5) can emphasize on precision more than recall, i.e. reduce FPR
  
  library(ROCR)
  library(survival)
  library(survcomp)
  
  ppv_measure <- function(label, pred.class, ytst){
    ifelse(all(pred.class!=label), 
           0, sum(pred.class==label & ytst==label)/sum(pred.class==label))
  }
  
  fpr_measure <- function(label, pred.class, ytst){
    sum(pred.class==label & ytst!=label)/sum(ytst!=label)
  }
  
  tpr_measure <- function(label, pred.class, ytst){
    sum(pred.class==label & ytst==label)/sum(ytst==label)
  }
  
  f_measure <- function(b, precision, recall){
    ifelse(precision == 0 && recall == 0,
           0, (1+b*b)*precision*recall/(b*b*precision+recall))
  }
  
  # error control
  if (is.null(pred$class))
    return(NULL)
  
  if(!is.factor(pred$class) && !is.character(pred$class)){
    stop("Indefinite predicted response type!")
  }
  
  if(length(pred$class) != length(ytst)){
    stop("Inconsistent prediction length!")
  }
  
  # remove NAs in ytst
  if(NA %in% ytst){
    id <- is.na(ytst)
    ytst <- ytst[!id]
    pred$class <- pred$class[!id]
    pred$prob <- pred$prob[!id, ]
  }
  
  # quantitative criteria
  acc <- mean(pred$class==ytst, na.rm=TRUE)
  ppv <- sapply(pos.label, ppv_measure, pred.class = pred$class, ytst = ytst)
  names(ppv) <- pos.label
  fpr <- sapply(pos.label, fpr_measure, pred.class = pred$class, ytst = ytst)
  names(fpr) <- pos.label
  tpr <- sapply(pos.label, tpr_measure, pred.class = pred$class, ytst = ytst)
  names(tpr) <- pos.label
  
  # f-score
  fval <- f_measure(b=beta,precision=ppv,recall=tpr)
  
  # concordance index for survival risk prediction
  if(!is.null(ysurv) && length(pos.label) == 1){
    ci <- concordance.index(pred$prob[ , pos.label, drop=TRUE], 
                            surv.time=ysurv[,1], surv.event=ysurv[,2])$c.index
  } else{
    ci <- NA
  }
  
  # auroc
  if(length(pos.label) == 1){
    ytst.bin <- as.character(ytst)
    ytst.bin[ytst == pos.label] <- paste0("pos")
    ytst.bin[ytst != pos.label] <- paste0("neg")
    ytst.bin <- factor(ytst.bin, levels = c("neg","pos"), ordered = TRUE)
    stopifnot(!any(is.na(ytst.bin)))
    predROCR <- prediction(pred$prob[ , pos.label, drop=TRUE], ytst.bin)
    auroc <- as.numeric(performance(predROCR, "auc")@y.values)
    if(length(auroc)>1){
      warning("Average AUC score computed!")
      auroc <- mean(auroc, na.rm=TRUE)
    }
    # plot ROC curve
    if(!is.null(savepath)){
      plotROCcv(res=list(test.prob=pred$prob, true.class=ytst.bin, pos.label=pos.label), 
                savepath=savepath, beta = beta)
    }
  } else{
    auroc <- NA
  }
  
  # return results
  return(list(acc=acc, fpr=fpr, tpr=tpr, ppv=ppv, fval=fval, concordance.index=ci, auroc=auroc))
}



indepValidation <- function(xtr, ytr, xtst, ytst, predictor, ysurv = NULL, 
                            pos.label = tail(names(table(ytst)),1), 
                            remove.const = TRUE, pthres = 0, 
                            ..., save.model = TRUE, seed = 53359292)
{
  # xtr, xtst are n*p feature matrices and ytr, ytst are label vectors
  # predictor is a char denoting the predictor function name (starting with '^predictor.*')
  # pos.label is a single or a vector of characters indicating the positive label(s) in ytst, defaulted the last level in ytst
  # remove.const defaulted TRUE to remove constant feature columns in xtr
  # pthres defaulted 0 does no indep signif feature reduction
  # indepValidation() returns model performances validated against test set which have been trained on training set
  
  library(survival)
  
  if (!is.null(seed)) set.seed(seed)
  
  # adjust correp rows and cols position in case they do not match
  if (!is.null(rownames(xtr)) && !is.null(rownames(xtst)) && !is.null(names(ytr)) && !is.null(names(ytst))) {
    samplelist.tr <- intersect(rownames(xtr), names(ytr))
    xtr <- xtr[samplelist.tr, ,drop=F]
    ytr <- ytr[samplelist.tr]
    
    samplelist.tst <- intersect(rownames(xtst), names(ytst))
    xtst <- xtst[samplelist.tst, ,drop=F]
    ytst <- ytst[samplelist.tst]
  }
  if (!is.null(colnames(xtr)) && !is.null(colnames(xtst))) {
    featlist <- intersect(colnames(xtr), colnames(xtst))
    xtr <- xtr[ ,featlist,drop=F]
    xtst <- xtst[ ,featlist,drop=F]
  } else {
    featlist <- colnames(xtr)
  }
  
  # remove NAs in ytr
  if(NA %in% ytr){
    id <- is.na(ytr)
    ytr <- ytr[!id]
    xtr <- xtr[!id, ]
  }
  
  # remove constant features
  if (remove.const) {
    d <- removeConst(xtr=xtr, xtst=xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # indep signif test for feature reduction
  if (pthres > 0) {
    plist <- featselectIndepSignif(xtr=xtr, ytr=ytr, pthres = pthres, pos.label = pos.label, ...)
    featlist <- names(plist)
    if (length(featlist) == 0) {
      stop("No features survived after independence significance test hence no feature reduction is reserved!")
    }
    xtr <- xtr[ ,featlist,drop=F]
    xtst <- xtst[ ,featlist,drop=F]
  }
  
  # model and predict and record system time as well
  
  pt <- proc.time()
  if(length(unique(as.character(ytr))) == 1){
    warning("Singe-class for training set!")
    pred <- predictorConstant(xtr, xtst, ytr, ...)
  } else{
    if(!is.character(predictor))
      stop("Specify predictor with a character!")
    FUN <- get(predictor, mode = "function")
    pred <- FUN(xtr=xtr, xtst=xtst, ytr=ytr, ...)
  }
  pt <- (proc.time() - pt)[["user.self"]]
  
  # evaluation and result list
  res <- c(list(validation="seqVal", predictor=predictor, cutoff=pred$cutoff, featlist=featlist, 
                true.surv=ysurv, true.class=ytst, test.class=pred$class, test.prob=pred$prob, 
                pos.label = pos.label, system_time = pt), 
           evaluatePred(pred=pred, ytst=ytst, ysurv=ysurv, pos.label = pos.label, ...))
  
  if(save.model){
    res$model <- pred$model
  }
  
  return(res)
}



crossValidation <- function(xtr, ytr, ..., seed=61215942, nfolds=5, nrepeats=10)
{
  # nfolds, nrepeats are used to generate nrepeats simultaneous training splits each being nfolds
  # crossValidation() returns cross validated results from a single dataset
  
  # data split
  if (!is.null(seed)) set.seed(seed)
  foldIndices <- caret::createMultiFolds(1:nrow(xtr), k=nfolds, times=nrepeats)
  
  message(nrepeats," repeated experiments of ",nfolds, "-fold cross validation!")
  foldres <- lapply(foldIndices, function(fold){
    message("New fold started...")
    ivres <- indepValidation(xtr=xtr[fold,,drop=F], ytr=ytr[fold], xtst=xtr[-fold,,drop=F], ytst=ytr[-fold], ..., save.model=FALSE)
    return(ivres)
  })
  
  cvres <- crossValidationCombineResults(foldres)
  
  return(cvres)
}



crossValidationCombineResults <- function(foldres, 
                                          elem_ave=c("acc","fpr","tpr","ppv","fval","concordance.index","auroc","system_time"), 
                                          elem_keepfold=c("true.surv","true.class","test.class","test.prob","featlist"),
                                          elem_inherit=c("predictor", "cutoff", "pos.label"))
{
  # be careful to use this function only with consistent cv foldres applied to the same xname,yname,prname,etc
  # elem_ave are those scores to average over cv fold
  # elem_keepfold are those to keep as they come from each fold
  # elem_inherit are those to inherit some parameters from res of subfold as in ensembleResults()
  
  cvres <- list("validation" = "cv_combined")
  
  if (is.null(foldres) || (is.list(foldres) && length(foldres) == 0)) {
    cvres <- list()
    for (elem in c(elem_inherit, elem_keepfold, elem_ave)) {
      cvres[[elem]] <- NA
    }
  } else {
    elem_red <- setdiff(names(foldres[[1]]), c(elem_inherit, elem_keepfold, elem_ave, "model", "validation"))
    if (length(elem_red) > 0) {
      warning("following elements are not recognized by default and are kept as in folds: ", 
              paste(elem_red, collapse = " , "))
      elem_keepfold <- setdiff(names(foldres[[1]]), c(elem_inherit, elem_ave))
    }
    cvres <- c(cvres,foldres[[1]][elem_inherit])
    for(elem in elem_keepfold){
      cvres[[elem]] <- lapply(foldres, '[[', elem)
    }
    for(elem in elem_ave){
      ss <- sapply(foldres, '[[', elem)
      if (elem == "ppv" || elem == "fval")
        ss[is.na(ss)] <- 0
      cvres[[elem]] <- mean(ss, na.rm=TRUE)
    }
  }
  
  return(cvres)
}



plotROCcv <- function(res, savepath, pos.label = tail(res$pos.label, 1), beta = 1)
{
  # res is the result list from calling indepValidation() or crossValidation(), ess a list having the following as elements 
  # true.class a vector of predicted class labels (can be cross-fold result given by list)
  # test.prob a matrix of predicted prob for each obs (in rows) indicating membership prob belonging to each class (in cols), or a list of such matrices for cv results
  # pos.label giving positive labels, must be a single character string for ROC to plot in this function
  
  library(ROCR)
  
  if (is.null(res) || (is.list(res) && length(res) == 0))
    return(NULL)
  
  # convert res$test.prob to a vector containing only pos.label
  stopifnot(length(pos.label) == 1)
  if (is.matrix(res$test.prob)) {
    res$test.prob <- res$test.prob[ , pos.label, drop = TRUE]
  } else if (is.list(res$test.prob)) {
    res$test.prob <- lapply(res$test.prob, function(u){
      stopifnot(is.matrix(u))
      u[ , pos.label, drop = TRUE]
    })
  } else {
    stop("res$test.prob should be a matrix or a list of matrix")
  }
  
  alpha <- function(b){return(1/(1+b*b))}
  
  if (!grepl("[.]pdf$", savepath)) savepath <- paste0(savepath, "_ROC_.pdf")
  
  pdf(savepath, height=10, width=15)
  
  predROCR <- prediction(res$test.prob, res$true.class)
  
  par(mfrow=c(2,2))
  
  acccurve <- performance(predROCR, "acc")
  plot(acccurve, lwd=2, main="Classification accuracy")
  grid()
  
  fcurve <- performance(predROCR, "f", alpha=alpha(b=beta))
  plot(fcurve, lwd=2, main=paste0("F-measure (beta=",beta,")"))
  grid()
  
  roccurve <- performance(predROCR, "tpr", "fpr")
  plot(roccurve, lwd=2, main="ROC curve")
  grid()
  
  # 	tprcurve <- performance(predROCR, "tpr")
  # 	plot(tprcurve, lwd=2, main="True Positive Rate (Recall or Sensitivity)")
  # 	grid()
  
  # 	fprcurve <- performance(predROCR, "fpr")
  # 	plot(fprcurve, lwd=2, main="False Positive Rate (Fallout)")
  # 	grid()
  
  # 	liftcurve <- performance(predROCR, "lift", "rpp")
  # 	plot(liftcurve, lwd=2, main="Lift Chart")
  # 	grid()
  
  # 	ppvcurve <- performance(predROCR, "ppv")
  # 	plot(ppvcurve, lwd=2, main="Positive Predictive Value (Precision)")
  # 	grid()
  
  prcurve <- performance(predROCR, "prec", "rec")
  plot(prcurve, lwd=2, main="Precision/Recall Graph")
  grid()
  
  dev.off()
  
  message("plot saved to ", dQuote(savepath))
  return(savepath)
}



# feature reduction -------------------------------------------------------

# in this section make sure named after featselectXXX where XXX comes from predictorXXX
# @param see predictorXXX
# @param keep.signif should only significant features be returned
# @return A vector named by corresponding features recording appropriate values sorted by decreasing "importance"



featselectIndepSignif <- function(xtr, ytr, pthres = 0.05, test = "t.test", method = "BH", keep.signif = TRUE,
                                  pos.label = tail(names(table(ytr)),1), ...)
{
  # xtr should be assigned col names and ytr is a vector of labels
  # test denotes a char from c("t.test","wilcox.test",...)
  # method denotes multiple test correction, default "BH", see ?p.adjust.methods
  # pos.label is a single character so that in case ytr has more than 2 levels we create a binary label wrt pos.label to perform the test
  # featselectIndepSignif() runs feature-by-feature independence significance test between contrasting classes
  # and returns p-values for each feature
  
  # remove const
  xtr <- removeConst(xtr = xtr, ...)
  # align samples
  if (!is.null(rownames(xtr)) && !is.null(names(ytr))) {
    samplelist <- intersect(rownames(xtr), names(ytr))
    xtr <- xtr[samplelist, ,drop=F]
    ytr <- ytr[samplelist]
  }
  # in case ytr has more than 2 levels
  stopifnot(length(pos.label) == 1)
  ytr.bin <- as.character(ytr)
  ytr.bin[ytr == pos.label] <- paste0("pos")
  ytr.bin[ytr != pos.label] <- paste0("neg")
  ytr.bin <- factor(ytr.bin, levels = c("neg","pos"), ordered = TRUE)
  
  # perform test
  featlist <- colnames(xtr)
  plist <- sapply(featlist, function(featname){
    fo <- paste0(featname," ~ y")
    t <- get(test)(formula = as.formula(fo), 
                   data = data.frame(y=as.factor(as.character(ytr.bin)),as.data.frame(xtr)), 
                   ...)
    t$p.value
  })
  names(plist) <- featlist
  plist <- p.adjust(p = plist, method = method)
  
  if (keep.signif)
    plist <- plist[plist < pthres]
  plist <- sort(plist, decreasing = FALSE)
  
  return(plist)
}



featselectRF <- function(model = NULL, ..., vthres = 0, keep.signif = TRUE)
{
  # computes importance measure as mean decrease in accuracy (type=1) for each feature
  # returns only those with positive importance ordered by decreasing value
  # NOTE model must be fit with randomForest() and with importance set TRUE
  
  library(randomForest)
  
  if (is.null(model))
    model <- predictorRF(...)$model
  s <- drop(randomForest::importance(x = model, type = 1, ...))
  
  if (keep.signif)
    s <- s[s > vthres]
  s <- sort(s, decreasing = TRUE)
  
  return(s)
}



featselectLogitLasso <- function(model = NULL, ..., keep.signif = TRUE)
{
  # computes absolute value of coefficients for each feature
  # returns only those with non-zero contribution to at least one of the phenotypes ordered by decreasing sum contribution
  # NOTE model must be fit with cv.glmnet() instead of glmnet() and only makes sense when do.normalize set TRUE as we do sum contribution
  
  library(glmnet)
  
  if (is.null(model))
    model <- predictorLogitLasso(...)$model
  s <- predict(object = model, newx = NULL, type = "coefficients", ...)
  s <- do.call('cbind', s)
  s <- rowSums(abs(s[-1, ])) # remove intercept
  
  if (keep.signif)
    s <- s[s > 0] # can be strongly positive or negative related!
  s <- sort(s, decreasing = TRUE)
  
  return(s)
}



featselectPAM <- function(model = NULL, ..., keep.signif = TRUE)
{
  # computes shrunken centroids in absolute values of coefficients for each feature
  # returns only those with non-zero contribution to at least one of the phenotypes ordered by decreasing sum contribution
  # NOTE model must be fit with pamr.train() and only makes sense when do.normalize set TRUE as we do sum contribution
  
  library(pamr)
  
  if (is.null(model))
    model <- predictorPAM(...)$model
  i <- which(model$threshold == model$best.thres)
  centroid <- pamr::pamr.predict(fit = model, newx = NULL, threshold = model$best.thres, type = "centroid", ...)
  
  # the following is adapted according to pamr::pamr.predict()
  delta.shrunk <- (centroid - model$centroid.overall)/model$sd
  s <- drop(abs(delta.shrunk) %*% rep(1, ncol(model$centroids)))
  
  if (keep.signif)
    s <- s[s > 0]
  s <- sort(s, decreasing = TRUE)
  
  return(s)
}



# classifiers -------------------------------------------------------------

# in this section make sure named after predictorXXX
# @param xtr/xtst Typically matrices that have ntr/ntst samples in rows, p explanatory variables in cols.
# @param ytr A ntr-vector of (binary or multi) class labels.
#   In implementation, "ytr" is always first converted to character in case that "ytr" has some factor level with zero obs.
# @param cutoff A numeric value only effective for binary classification such that
#   - except for SparseSVM, "cutoff" is a threshold between 0 and 1, defaulted 0.5, to cut predicted prob;
#   - for SparseSVM, "cutoff" is a threshold between -Inf and Inf, defaulted 0, to cut margin decision values.
# @param ... Other algorithm-specific arguments to be passed onto the main function within.
# @return A list with entries:
# \item{model}{The trained model with "xtr" and "ytr"}
# \item{class}{A character vector of predicted classes}
# \item{prob}{A numeric matrix of dimension nrow(xtst) X nlevels(ytr) such that
#   - except for SparseSVM, "prob" records predicted prob (may contain NAs for current implementation for multi-class KNN)
#   - for SparseSVM, "prob" records predicted prob margin decision values}
# \item{cutoff}{Simply a copy of the input parameter "cutoff"}



predictorLDA <- function(xtr, xtst, ytr, cutoff = 0.5, do.normalize = TRUE, ...)
{
  # no feature selection, no parameter tuning
  
  library(MASS)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # train
  model <- lda(x = xtr, grouping = as.factor(as.character(ytr)), ...)
  # predict
  pred <- predict(object = model, newdata = xtst, ...)
  pred.class <- as.character(pred$class)
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred$posterior)] <- pred$posterior
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorLogitLasso <- function(xtr, xtst, ytr, alpha = 1, cutoff = 0.5, do.normalize = TRUE, ...)
{
  # Logit regression with elastic net penalty, default alpha=1, which should not be altered, denotes lasso for feature selection
  # parameter tuning as implemented by cv.glmnet()
  
  library(glmnet)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  # train
  model <- cv.glmnet(x = as.matrix(xtr), y = as.factor(as.character(ytr)), family = "multinomial", 
                     alpha = alpha, standardize = do.normalize, ...)
  # predict
  pred <- drop(predict(object = model, newx = as.matrix(xtst), type = "response", ...))
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorLinearSVM <- function(xtr, xtst, ytr, 
                               kernel = "linear", do.normalize = TRUE, 
                               cost = 10^(-3:3), cross = 5, cutoff = 0.5, ...)
{
  # cost is a vector of C parameter grid to tune with, no feature selection
  # NOTE predictorLinearSVM, though default to linear kernel which should not be altered, is a wrapper function for any kernel methodsm, see predictorRadialSVM()
  # NOTE cross validation is set to always minimize misclassification error by default of tune.svm arguments despite that it might not be appropriate for largely unbalanced classes
  
  library(e1071)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  # train
  model <- tune.svm(x = xtr, y = as.factor(as.character(ytr)), kernel = kernel, 
                    scale = do.normalize, type = "C-classification", probability = TRUE, 
                    cost = cost, tunecontrol = tune.control(sampling="cross", cross=cross), ...)
  # predict
  pred <- predict(object = model$best.model, newdata = xtst, probability = TRUE, ...)
  pred.class <- as.character(pred)
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(attr(pred,"probabilities"))] <- attr(pred,"probabilities")
  
  res <- list(model=model$best.model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorRadialSVM <- function(xtr, xtst, ytr, do.normalize = TRUE, ...)
{
  # implemented by replacing linear kernel by RBF kernel in predictorLinearSVM()
  # gamma for radial kernel is set with median trick implemented by kernlab:sigest()
  
  gamma <- kernlab::sigest(xtr, scaled = do.normalize)['50%']
  res <- predictorLinearSVM(xtr = xtr, xtst = xtst, ytr = ytr, kernel = "radial", gamma = gamma, do.normalize = do.normalize, ...)
  return(res)
}



predictorKendallSVM <- function(xtr, xtst, ytr, kernel = pcaPP::cor.fk, do.normalize = FALSE, 
                                kmat = NULL, 
                                cost = 10^(-3:3), cross = 5, cutoff = 0.5, ...)
{
  # by the concept with default kendall kernel, do.normalize should be set FALSE by default as well!! Be careful when calling with other kernels
  # kernel is provided as a "function" instead of a "character"
  # kmat is used if provided to accelerate the computation and names must be matched to rows in xtr/xtst to identity tr/tst indices
  # cost is a vector of C parameter grid to tune with, no feature selection
  # NOTE cross validation is set to always minimize misclassification error by default of tune.svm arguments despite that it might not be appropriate for largely unbalanced classes
  # NOTE this function also is a generic function which uses kernlab implementation and can thus deal with explicit kernel matrix as input with any kernel function 
  #     and it can be slower for certain kernels (and how cross is done) compared to e1071 implementation (linear kernel for instance)
  
  library(kernlab)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if (is.null(kmat)) {
    if(do.normalize){
      d <- normalizeData(xtr, xtst, ...)
      xtr <- d$xtr
      xtst <- d$xtst
    }
    stopifnot(is.function(kernel))
    kmat <- computeKernelMatrix(x = rbind(xtr, xtst), kernel = kernel, ...)
  } else {
    stopifnot(!is.null(xtr) && !is.null(xtst) && !is.null(dimnames(kmat)))
    idx <- match(c(rownames(xtr),rownames(xtst)), rownames(kmat))
    kmat <- kmat[idx,idx]
  }
  
  # CV for tuning cost
  train.idx <- 1:nrow(xtr)
  cost <- sort(cost, decreasing = FALSE) # parsimony principle
  cv.acc <- sapply(cost, function(cpm){
    cv.model <- ksvm(x = as.kernelMatrix(kmat[train.idx, train.idx]), 
                     y = as.factor(as.character(ytr)), 
                     scaled = do.normalize, kernel = kernel, type = "C-svc", 
                     C = cpm, cross = cross, prob.model = FALSE, ...)
    cv.model@cross
  })
  best.cpm <- cost[which.min(cv.acc)]
  # train
  model <- ksvm(x = as.kernelMatrix(kmat[train.idx, train.idx]), 
                y = as.factor(as.character(ytr)), 
                scaled = do.normalize, kernel = kernel, type = "C-svc", 
                C = best.cpm, cross = 0, prob.model = TRUE, ...)
  # predict
  pred <- predict(object = model, 
                  newdata = as.kernelMatrix(kmat[-train.idx,train.idx,drop=F][,SVindex(model),drop=F]), 
                  type = "probabilities", ...)
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorKNN <- function(xtr, xtst, ytr, k = seq(1,10,2), do.normalize = TRUE, cutoff = 0.5, ...)
{
  # no feature selection, with parameter tuning for k
  # NOTE k takes only odd values to avoid confusion from ties
  # NOTE cross validation is set to minimize misclassification error by default of tune.knn arguments
  # IMPORTANT for multi-class (>2) case the prob of winning class is the proportion of votes of winning class from neighbours, 
  #           and the other classes can take NA values as class::knn has not yet implemented this option
  
  library(class)
  library(e1071)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # train
  model <- tune.knn(x = xtr, y = as.factor(as.character(ytr)), k = k, 
                    tunecontrol=tune.control(sampling="cross", cross=5), ...)
  best.k <- model$best.parameters$k
  # predict
  pred <- knn(train = xtr, test = xtst, cl = ytr, k = best.k, prob = TRUE, ...)
  pred.class <- as.character(pred)
  names(pred.class) <- rownames(xtst)
  # IMPORTANT for multi-class (>2) case NA denotes unknown values!!
  pred.prob <- matrix(ifelse(length(classes.eff) == 2, 0, NA), 
                      nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[attr(pred,"prob") == 1, ] <- 0
  id.col <- match(pred.class, classes)
  pred.prob[cbind(1:nrow(pred.prob), id.col)] <- attr(pred,"prob")
  # in case of two effective classes (with at least one obs in such class), fill out the losing.prob by 1-winning.prob
  if (length(classes.eff) == 2) {
    pred.prob[cbind(1:nrow(pred.prob), sum(match(classes.eff, classes)) - id.col)] <- 1 - attr(pred,"prob")
  }
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorNB <- function(xtr, xtst, ytr, do.normalize = TRUE, cutoff = 0.5, ...)
{
  # no feature selection, no parameter tuning
  
  library(e1071)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # train
  model <- naiveBayes(x = xtr, y = as.factor(as.character(ytr)), ...)
  # predict
  pred <- predict(object = model, newdata = xtst, type = "raw", ...)
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorGBM <- function(xtr, xtst, ytr, n.trees = 1500, shrinkage = 0.002, 
                         interaction.depth = 6, bag.fraction = 1, cutoff = 0.5, ...)
{
  # no feature selection
  # NOTE no need for do.normalize for base learner DT
  # NOTE tuning for algorithmic parameters is not yet implemented
  
  library(gbm)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  # train
  model <- gbm.fit(x = xtr, y = as.factor(as.character(ytr)), 
                   n.trees = n.trees, shrinkage = shrinkage, interaction.depth = interaction.depth, 
                   bag.fraction = bag.fraction, distribution = "multinomial", ...)
  # predict
  pred <- drop(predict(object = model, newdata = xtst, n.trees = model$n.trees, type = "response", ...))
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorRF <- function(xtr, xtst, ytr, ntree = 1500, importance = TRUE, cutoff = 0.5, ...)
{
  # ntree is a vector of parameter grid to tune with, random feature selection
  # NOTE cross validation is set to always minimize misclassification error by default of tune.svm arguments despite that it might not be appropriate for largely unbalanced classes
  # NOTE no need for do.normalize for base learner DT
  
  library(randomForest)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  # train
  model <- randomForest(x = xtr, y = as.factor(as.character(ytr)), 
                        ntree = ntree, importance = importance, ...)
  # predict
  pred <- predict(object = model, newdata = xtst, type = "prob", ...)
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorSparseSVM <- function(xtr, xtst, ytr, cost = 1, cutoff = 0, do.normalize = TRUE, ...)
{
  # cost is a single C parameter without being tuned
  # LiblineaR(type=5) L1-regularized L2-loss support vector classification (i.e. with feature selection)
  # NOTE tuning for algorithmic parameters is not yet implemented
  # IMPORTANT currently the returned entry "prob" in the result list does NOT mean prob indeed as they are in fact margin decision values,
  #           and thus "cutoff" takes 0 by default and should range over -Inf to Inf,
  #           as LiblineaR has not yet implemented this option but practically doable, in reference to
  #           Wu, T.F., Lin, C.J. and Weng, R.C., "Probability estimates for multi-class classification by pairwise coupling." JMLR 5 (2004): 975-1005.
  
  library(LiblineaR)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # train
  # in case of two-class case, ensure that straight zeros are stored for the first case while
  # non-zero margin values are recorded for the second class (considered larger or positive by default)
  o <- order(ytr, decreasing = TRUE)
  model <- LiblineaR(data = xtr[o, , drop = FALSE], target = as.factor(as.character(ytr))[o], 
                     type = 5, cost = cost, ...)
  # predict
  pred <- predict(object = model, newx = xtst, decisionValues = TRUE, ...)
  pred.class <- as.character(pred$predictions)
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred$decisionValues)] <- pred$decisionValues
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorPAM <- function(xtr, xtst, ytr, do.normalize = TRUE, cross = 5, cutoff = 0.5, ...)
{
  # known as nearest shrunken centroid classifier
  # feature selection by thresholding centroids
  # parameter tuning for best threshold implemented by pamr.cv()
  # NOTE returned model has one more entry recording best threshold from cv runs
  
  library(pamr)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  # train
  dat <- list(x = t(xtr), y = as.factor(as.character(ytr)))
  model <- pamr.train(data = dat, ...)
  model.cv <- pamr.cv(fit = model, data = dat, nfold = cross, ...)
  model$best.thres <- model.cv$threshold[which.max(model.cv$loglik)]
  # predict
  pred <- pamr.predict(fit = model, newx = t(xtst), 
                       threshold = model$best.thres, type = "posterior", ...)
  pred.class <- colnames(pred)[max.col(pred)]
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  pred.prob[ , colnames(pred)] <- pred
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



predictorConstant <- function(xtr, xtst, ytr, cutoff = 0.5, ...)
{
  # predict ytst by a constant class as the modal class in ytr and a constant frequency as the proportion present in ytr
  # NOTE indeed a baseline predictor and an error control for other predictors when only one class is present effectively
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  # stopifnot(length(classes.eff) >= 2) # also work with single class in training data
  modal <- classes[which.max(tab)]
  
  model <- "Modal-class classifier"
  pred.class <- rep(modal, nrow(xtst))
  names(pred.class) <- rownames(xtst)
  pred.prob <- matrix(tab/length(ytr), byrow = TRUE, 
                      nrow = nrow(xtst), ncol = length(classes), 
                      dimnames = list(rownames(xtst), classes))
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}



# utiles ------------------------------------------------------------------



# NA proportion function
NArate <- function(v)
{
  mean(is.na(v))
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#   - cols:   Number of columns in layout
#   - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL)
{
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Taken from http://www.r-bloggers.com/memory-management-in-r-a-few-tips-and-tricks/
# improved list of objects

.ls.objects <- function(pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5)
{
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10)
{
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


