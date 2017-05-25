# Randomly select features, for later test if the FS truely improves classification accuracy
# The classifiers below will be tested: SVM, RIDGE, LASSO, Random Forest
# Evaluation will be based on AUC (with probability)
# Selected features will not be recorded, but stability will be calculated

remove(list = ls())
time.start <- Sys.time()
print(paste('Start time:', time.start))

library(xlsx)
library(randomForest)
library(class)
library(e1071)
library(glmnet)
library(pROC)
library(reshape)
library(ggplot2)
library(cowplot)
library(infotheo)

stab <- function(List){
  nsample <- length(List)  # times of resampling, should be the number of samples, if the main function uses LOOCV
  num.select.genes <- length(List[[1]])
  union <- Reduce(union, List)  # F in the formula: the list of all features, which have been selected in at least one of n sampling steps
  whole <- unlist(List)  # all the features in List (with replicate)
  sum <- 0
  for(fea in union){
    freq <- length(which(whole == fea))
    sum <- sum + freq
  }
  stab.na <- sum / (nsample * length(union))
  return(stab.na)
}

# prepare for the data set
data.raw <- read.xlsx("mercury.xlsx", sheetIndex = 2, colNames = TRUE)  # each row is a protein, each column is a sample, but the first two columns are pretein names and id

nsample <- 19  # number of samples
nfeature <- nrow(data.raw) - 1  # number of features (proteins)
num.control <- 10
num.high <- 9
num.low <- 9
# num.medium <- 10

x.raw <- log2(data.raw[ , c(3:12, 13:21)] + 1)   # nsample samples (control group and high dose group)
y.raw <- c(rep(1,num.control), rep(2,num.high))  # control samples are labeled as 1, and treatment samples are labeled as 2
pr.name <- data.raw$accession
pr.id <- data.raw$id

step <- 10
num.select.genes.down <- 10
num.select.genes.up <- 400
# num.select.genes <- 100
num.resample <- 1000
# num.auc <- num.resample

auc.rdm.rf.list <- list()
auc.rdm.svm.list <- list()
auc.rdm.ridge.list <- list()
auc.rdm.lasso.list <- list()

stab.rdm.list <- list()

# loop of number of selected genes
count.loop.sel <- 0
for (num.select.genes in seq(num.select.genes.down, num.select.genes.up, step)){
  print(paste('num of selected genes', num.select.genes))
  count.loop.sel <- count.loop.sel + 1
  
  auc.rdm.rf.list[[count.loop.sel]] <- list()
  auc.rdm.svm.list[[count.loop.sel]] <- list()
  auc.rdm.ridge.list[[count.loop.sel]] <- list()
  auc.rdm.lasso.list[[count.loop.sel]] <- list()
  
  stab.rdm.list[[count.loop.sel]] <- list()
  
  seed <- 100  # seed for resampling
  # loop of resampling
  for (count.loop in seq(1, num.resample)){
    seed <- seed + 1  # change the seed every time, to avoid same resampling
    
    select.genes.rdm.list <- list()
    
    pred.rdm.rf.list <- rep(0, nsample)  # predictions of all samples in this loop of feature #
    
    pred.rdm.svm.list <- rep(0, nsample)
    
    pred.rdm.ridge.list <- rep(0, nsample)
    
    pred.rdm.lasso.list <- rep(0, nsample)
    
    # LOOV
    for (i in c(1:nsample)){
      seed <- seed + 1
      x.train.temp <- x.raw[, -i]
      y.train <- y.raw[-i]
      x.test.temp <- x.raw[, i]
      y.test <- y.raw[i]
      
      # Normalisation
      x.mean <- apply(x.train.temp, 1, mean)
      x.sd <- apply(x.train.temp, 1, sd)
      x.train <- (x.train.temp - x.mean) / x.sd  # row --> feature, column --> sample
      x.test <- (x.test.temp - x.mean) / x.sd
      
      # randomly select genes
      id.genes <- seq(1, nfeature)
      set.seed(seed)
      select.genes.rdm <- sample(id.genes, num.select.genes, replace = FALSE)
      select.genes.rdm.list[[i]] <- select.genes.rdm
      
      # Data sets for classifiers, usually row --> sample, column --> feature
      
      train.rdm <- t(x.train[select.genes.rdm, ])
      label.train <- y.train
      
      test.rdm <- t(x.test[select.genes.rdm])
      colnames(test.rdm) <- colnames(train.rdm)
      label.test <- y.test
      
      # use Random Forest to test the performance of three feature subsets
      label.train.rf <- factor(label.train)
      label.test.rf <- factor(label.test)
      
      set.seed(1234)
      model.rdm.rf <- randomForest(train.rdm, label.train.rf)
      
      pred.rdm.rf <- predict(model.rdm.rf, test.rdm, type = "prob")[2]
      
      pred.rdm.rf.list[i] <- pred.rdm.rf
      
      # use SVM to test the performance of three feature subsets
      train.rdm.svm <- as.data.frame(cbind(train.rdm, label.train), col.names = c(colnames(train.rdm), 'label.train'))
      test.rdm.svm <- as.data.frame(test.rdm)
      
      model.rdm.svm <- svm(formula = factor(label.train)~., data = train.rdm.svm, probability = TRUE)
      
      pred.rdm.svm <- predict(model.rdm.svm, test.rdm.svm, probability = TRUE)
      
      pred.rdm.svm.list[i] <- attr(pred.rdm.svm, "probabilities")[2]
      
      # RIDGE REGRESSION 
      ridge.cv.rdm <- cv.glmnet(train.rdm, label.train, alpha = 0, type.measure="mse",
                                family = "binomial", standardize=F, nfolds = 3)
      # Optimal lambda:
      ridge.lam.rdm <- ridge.cv.rdm$lambda.min
      
      # Determine predicted values
      pred.rdm.ridge <- predict(ridge.cv.rdm, s = ridge.lam.rdm, type="response", newx = test.rdm)
      
      pred.rdm.ridge.list[i] <- pred.rdm.ridge
      
      # LASSO REGRESSION 
      lasso.cv.rdm <- cv.glmnet(train.rdm, label.train, alpha = 1, type.measure="mse",
                                family = "binomial", standardize=F, nfolds = 3)
      # Optimal lambda:
      lasso.lam.rdm <- lasso.cv.rdm$lambda.min
      
      # Determine predicted values
      pred.rdm.lasso <- predict(lasso.cv.rdm, s = lasso.lam.rdm, type="response", newx = test.rdm)
      
      pred.rdm.lasso.list[i] <- pred.rdm.lasso
    }
    
    
    # sumup the results
    auc.rdm.rf <- auc(y.raw, pred.rdm.rf.list)
    auc.rdm.svm <- auc(y.raw, pred.rdm.svm.list)
    auc.rdm.ridge <- auc(y.raw, pred.rdm.ridge.list)
    auc.rdm.lasso <- auc(y.raw, pred.rdm.lasso.list)
    
    auc.rdm.rf.list[[count.loop.sel]] <- append(auc.rdm.rf.list[[count.loop.sel]], auc.rdm.rf)
    auc.rdm.svm.list[[count.loop.sel]] <- append(auc.rdm.svm.list[[count.loop.sel]], auc.rdm.svm)
    auc.rdm.ridge.list[[count.loop.sel]] <- append(auc.rdm.ridge.list[[count.loop.sel]], auc.rdm.ridge)
    auc.rdm.lasso.list[[count.loop.sel]] <- append(auc.rdm.lasso.list[[count.loop.sel]], auc.rdm.lasso)
    
    stab.rdm <- stab(select.genes.rdm.list)
    stab.rdm.list[[count.loop.sel]] <- append(stab.rdm.list[[count.loop.sel]], stab.rdm)
  }
}

time.end <- Sys.time()
print(paste('End time:', time.end))
# # Fin