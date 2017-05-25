# This belongs to the GitHub project FScompare: https://github.com/zhxiaokang/FScompare 
# Study the effect of using different numbers of selected genes
# Use LOOCV to compare the performance of SAM, mRMR, and GeoDE, and use randomly selected genes as baseline
# The classifiers below will be tested: SVM, RIDGE, LASSO, Random Forest
# Evaluation will be based on AUC (with probability)
# Selected features will be recorded

remove(list = ls())

library(xlsx)
library(samr)
library(mRMRe)
library(randomForest)
library(class)
library(e1071)
library(glmnet)
library(pROC)
library(reshape)
library(ggplot2)
library(cowplot)
library(infotheo)

#  load package "GeoDE" but stop the in-function plotting
library(GeoDE)
chdirSig <- function(data,sampleclass,gammas=list(1.0),nnull=10,CalculateSig=FALSE)
{
  
  #
  # Calculate the pca
  #
  
  pca1 <- prcomp(t(as.matrix(data[-1])))
  
  meanvec <- rowMeans(as.matrix(data[-1][sampleclass==2]))-rowMeans(as.matrix(data[-1][sampleclass==1]))
  
  n1 <- sum(sampleclass==1)
  n2 <- sum(sampleclass==2)
  
  #
  #   In order to prevent singularities later we now select the number 
  #of PCs to keep
  #
  
  cumsum <- pca1$sdev^2/sum(pca1$sdev^2)
  keepPC <- length(cumsum[cumsum>0.001])
  
  #
  # Now calculte the characteristic direction in pca space
  #
  V <- pca1$rotation[,1:keepPC]
  R <- pca1$x[,1:keepPC]
  
  Dd <- (t(R[sampleclass==1,])%*%R[sampleclass==1,]+t(R[sampleclass==2,])%*%R[sampleclass==2,])/(n1+n2-2)
  
  sigma <- mean(diag(Dd))
  
  ShrunkMats <- lapply(gammas, function(x) solve(x*Dd + sigma*(1-x)*diag(keepPC)))
  
  # Map back to the full expression space
  
  
  b <- lapply(ShrunkMats, function(x) matrix(V%*%x%*%t(V)%*%meanvec,dimnames=list(c(as.list(as.character(data[[1]]))), 1)))
  
  # Normalize the characteristic directions
  
  
  b <- lapply(b, function(x) x/sqrt(sum(x^2)))
  
  # b<-as.vector(b[[1]])  
  #The projection of the characteristic direction in the first two PCs
  
  
  b2dscale <- colMeans(R[sampleclass==2,1:2])- colMeans(R[sampleclass==1,1:2])
  
  
  b2dscale <- sqrt(sum(b2dscale^2))
  
  #print(b2dscale)
  
  projchdir2d <-lapply(b, function(x) list( b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,1])), b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,2]))))  
  
  
  
  if (CalculateSig) 
  {
    
    
    ########################################
    # Generate a null distribution of chdirs
    ########################################
    
    
    
    
    nu<-n1+n2-2
    y1 <- t(t(mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd) *sqrt(nu / rchisq(nnull, nu))))
    y2 <- t(t(mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd) *sqrt(nu / rchisq(nnull, nu))))
    bmeanvec <- colMeans(R[sampleclass==2,])- colMeans(R[sampleclass==1,])
    y <- t(y2-y1)
    
    
    #  y <- mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd)
    #  y2<-sqrt(nu / rchisq(nnull, nu))
    
    ###############################################
    # For each value of gamma and each of the null
    ###############################################  
    #    print("1/2")
    pb <- txtProgressBar(min = 0, max = nnull,style=3) 
    nullchdirs <- lapply(gammas, function(x)
    {
      rowMeans( 
        as.data.frame(
          mapply( function(mv,count)
          {
            
            
            setTxtProgressBar(pb, count)
            
            
            sm <- solve(x*Dd + sigma*(1-x)*diag(keepPC))        
            bn <-  matrix(V%*%sm%*%t(V)%*%as.numeric(mv%*%t(V)),dimnames=list(c(as.list(as.character(data[[1]]))), 1))
            bn <- bn/sqrt(sum(bn^2))
            bn<-bn^2
            bn<-sort(bn,decreasing=TRUE)
            
            
          },as.data.frame(y),c(1:length(y))
          )
        )
      )
    }
    )
    close(pb)  
    
    
    #  print(class(nullchdirs))
    #  print(length(nullchdirs))
    #  print(dim(nullchdirs))
    
    ###########################################
    # The ratio to null
    ###########################################
    #    print("2/2")
    #    pb <- txtProgressBar(min = 0, max = nnull,style=3)
    
    
    ratio <- mapply(function(ba,bn,count) {
      #      setTxtProgressBar(pb, count)
      relerr <- sort(ba^2,decreasing=TRUE)/bn
      
      relerr <- cumsum(relerr)/sum(relerr)-c(1:length(meanvec))/length(meanvec)
    }, b,nullchdirs, c(1:length(nullchdirs)),SIMPLIFY=FALSE)
    #    close(pb)
    #    print("This may take several minutes, please wait...")
    nsiggenes <-lapply(ratio,function(x) which.max(x))
    
    list(chdir=b,pca2d=R[,1:2],chdir_pca2d=projchdir2d,null_rank_dist=nullchdirs,ratio_to_null=ratio,number_sig_genes=nsiggenes)
    
    
  }else
  {
    list(chdir=b,pca2d=R[,1:2],chdir_pca2d=projchdir2d) 
  } 
  
  
}

chdirAnalysis <- function(datain, sampleclass, gammas=list(1.0), nnull=10,CalculateSig=FALSE)
{  # Test that the inputs are in the correct form and are self-consisent
  
  if(length(sampleclass)!=(length(datain)-1)) stop("number of elements in sampleclass is inconsistent with input data")
  if(!is.data.frame(datain)) stop("Input data is not in the form of a data frame")
  if(FALSE%in%(c("1","2")%in%levels(sampleclass))) stop ("sample class does not include \'1\' and \'2\'")
  if(length(datain[sampleclass==1])<2) stop ("too few controll samples")
  if(length(datain[sampleclass==2])<2) stop ("too few samples")
  
  # Calculate the characteristic direction
  chdirresults <- chdirSig(datain,sampleclass,gammas,nnull=nnull,CalculateSig=CalculateSig)
  
  # produce plots
  #  print("plotting...")
  # chdirplots(chdirresults,sampleclass,gammas,CalculateSig)
  #  print("generating output...")
  
  # Generate result dataframe
  outAll <- lapply(chdirresults[[1]], function(x) {x[sort.list(x^2,decreasing=TRUE),]})
  if(CalculateSig)
  {
    outSig <- mapply( function(x,ns) {
      
      x[sort.list(x^2,decreasing=TRUE)[1:ns],]}, chdirresults[[1]],chdirresults[[6]],SIMPLIFY=FALSE)
    list(chdirprops=chdirresults,results=outSig)
  }else{list(chdirprops=chdirresults,results=outAll)}
  
}

# define the function of calculating stability
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
data.raw <- read.xlsx("TPM_Normalised_Counts_ALL.xlsx", "TPM")  # each row is a protein, each column is a sample, but the first two columns are pretein names and id

nsample <- 14  # number of samples
nfeature <- nrow(data.raw) - 1  # number of features (proteins)

x.raw <- log2(data.raw[ , c(3:10, 44:49)] + 1)   # nsample samples (control group and high dose group)
num.control <- 8
num.treat <- 6
if (nsample != num.control + num.treat){
  print('Number of samples do not agree!')
}
y.raw <- c(rep(1,num.control), rep(2,num.treat))  # control samples are labeled as 0, and treatment samples are labeled as 1
pr.name <- data.raw$accession
pr.id <- data.raw$id

step <- 2
num.select.genes.down <- 6
num.select.genes.up <- 10
num.auc <- ceiling((num.select.genes.up - num.select.genes.down + 1) / step)  # number of loops for different numbers of selected features

select.genes.sam.list <- list()
select.genes.mrmr.list <- list()
select.genes.geode.list <- list()
select.genes.rdm.list <- list()

stab.sam.list <- rep(0, num.auc)
stab.mrmr.list <- rep(0, num.auc)
stab.geode.list <- rep(0, num.auc)
stab.rdm.list <- rep(0, num.auc)

auc.sam.rf.list <- rep(0, num.auc)
auc.mrmr.rf.list <- rep(0, num.auc)
auc.geode.rf.list <- rep(0, num.auc)
auc.rdm.rf.list <- rep(0, num.auc)

auc.sam.svm.list <- rep(0, num.auc)
auc.mrmr.svm.list <- rep(0, num.auc)
auc.geode.svm.list <- rep(0, num.auc)
auc.rdm.svm.list <- rep(0, num.auc)

auc.sam.ridge.list <- rep(0, num.auc)
auc.mrmr.ridge.list <- rep(0, num.auc)
auc.geode.ridge.list <- rep(0, num.auc)
auc.rdm.ridge.list <- rep(0, num.auc)

auc.sam.lasso.list <- rep(0, num.auc)
auc.mrmr.lasso.list <- rep(0, num.auc)
auc.geode.lasso.list <- rep(0, num.auc)
auc.rdm.lasso.list <- rep(0, num.auc)

count.loop <- 0
# loop of number of selected genes
for (n in seq(num.select.genes.down, num.select.genes.up, step)){
  print(paste('num of selected genes', n))
  count.loop <- count.loop + 1
  num.select.genes <- n
  
  select.genes.sam.list[[count.loop]] <- list()
  select.genes.mrmr.list[[count.loop]] <- list()
  select.genes.geode.list[[count.loop]] <- list()
  select.genes.rdm.list[[count.loop]] <- list()
  
  pred.sam.rf.list <- rep(0, nsample)  # predictions of all samples in this loop of feature #
  pred.mrmr.rf.list <- rep(0, nsample)
  pred.geode.rf.list <- rep(0, nsample)
  pred.rdm.rf.list <- rep(0, nsample)
  
  pred.sam.svm.list <- rep(0, nsample)
  pred.mrmr.svm.list <- rep(0, nsample)
  pred.geode.svm.list <- rep(0, nsample)
  pred.rdm.svm.list <- rep(0, nsample)
  
  pred.sam.ridge.list <- rep(0, nsample)
  pred.mrmr.ridge.list <- rep(0, nsample)
  pred.geode.ridge.list <- rep(0, nsample)
  pred.rdm.ridge.list <- rep(0, nsample)
  
  pred.sam.lasso.list <- rep(0, nsample)
  pred.mrmr.lasso.list <- rep(0, nsample)
  pred.geode.lasso.list <- rep(0, nsample)
  pred.rdm.lasso.list <- rep(0, nsample)
  
  # LOOCV
  seed <- 1234  # a seed for "sample", need to change every loop of LOOCV, otherwise the "randomly" selected genes will all be the same
  for (i in c(1:nsample)){
    seed <- seed + 1
    
    x.train.temp <- x.raw[, -i]
    y.train <- y.raw[-i]
    x.test.temp <- x.raw[, i]
    y.test <- y.raw[i]
    
    # Normalisation
    #    x.mean <- apply(x.train.temp, 1, mean)
    #    x.sd <- apply(x.train.temp, 1, sd)
    #    x.train <- (x.train.temp - x.mean) / x.sd  # row --> feature, column --> sample
    #    x.test <- (x.test.temp - x.mean) / x.sd
    
    # Without normalisation
    x.train <- x.train.temp
    x.test <- x.test.temp
    
    # use SAM to select the significantly differentially regulated genes (proteins)
    data <- list(x=x.train,y=y.train,genenames=paste("gene",as.character(pr.id),sep=""), geneid=as.character(pr.id), logged2=TRUE)  
    invisible(capture.output(samr.obj <- samr(data, resp.type="Two class unpaired", nperms=100)))
    invisible(capture.output(delta.table <- samr.compute.delta.table(samr.obj)))
    del <- -3
    siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data, delta.table)
    siggenes.table$genes.lo[,4] <- -as.numeric((siggenes.table$genes.lo[,4]))  # reverse the Score(d) of genes.lo
    siggenes <- rbind(siggenes.table$genes.up, siggenes.table$genes.lo)  # combine the up and low genes
    siggenes <- siggenes[order(siggenes[,4],decreasing = TRUE),]  # sort the genes aucording to the Score(d)
    select.genes.sam <- siggenes[1:num.select.genes, 1]  # significantly differentially regulated genes selected by SAM
    select.genes.sam <- as.numeric(select.genes.sam) - 1
    select.genes.sam.list[[count.loop]] <- append(select.genes.sam.list[[count.loop]], list(select.genes.sam))
    
    # use Feature Selection method mRMR to select the important genes (proteins)
    data.mrmr <- cbind(y.train, t(x.train))  # in mRMRe, each row is a sample, and each colunm is a feature, and the first colunm is class label
    dd <- mRMR.data(data = data.frame(data.mrmr))
    result.mrmr <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 1, feature_count = num.select.genes)
    select.genes.mrmr <- result.mrmr@filters
    select.genes.mrmr <- as.numeric(select.genes.mrmr$`1`) - 1
    select.genes.mrmr.list[[count.loop]] <- append(select.genes.mrmr.list[[count.loop]], list(select.genes.mrmr))
    
    # use GeoDE to select the important genes (proteins)
    select.genes.geode <- rep(0, num.select.genes)
    gammas <- 1
    data.geode <- cbind(pr.name, x.train)
    chdir_analysis <- chdirAnalysis(data.geode, factor(y.train), gammas, CalculateSig=FALSE, nnull=10)
    select.genes.geode.names <- names(chdir_analysis$results[[1]][1:num.select.genes])
    for (s in seq(1, num.select.genes)){
      select.genes.geode[s] <- which(data.raw$accession == select.genes.geode.names[s])
    }
    select.genes.geode.list[[count.loop]] <- append(select.genes.geode.list[[count.loop]], list(select.genes.geode))
    
    # randomly select genes
    id.genes <- seq(1, nfeature)
    set.seed(seed)
    select.genes.rdm <- sample(id.genes, num.select.genes, replace = FALSE)
    select.genes.rdm.list[[count.loop]] <- append(select.genes.rdm.list[[count.loop]], list(select.genes.rdm))

    # Data sets for classifiers, usually row --> sample, column --> feature
    
    train.sam <- t(x.train[select.genes.sam, ])
    train.mrmr <- t(x.train[select.genes.mrmr, ])
    train.geode <- t(x.train[select.genes.geode, ])
    train.rdm <- t(x.train[select.genes.rdm, ])
    label.train <- y.train
    
    test.sam <- t(x.test[select.genes.sam])
    test.mrmr <- t(x.test[select.genes.mrmr])
    test.geode <- t(x.test[select.genes.geode])
    test.rdm <- t(x.test[select.genes.rdm])
    label.test <- y.test
    
    # use Random Forest to test the performance of three feature subsets
    label.train.rf <- factor(label.train)
    label.test.rf <- factor(label.test)
    
    set.seed(1234)
    model.sam.rf <- randomForest(train.sam, label.train.rf)
    set.seed(1234)
    model.mrmr.rf <- randomForest(train.mrmr, label.train.rf)
    set.seed(1234)
    model.geode.rf <- randomForest(train.geode, label.train.rf)
    set.seed(1234)
    model.rdm.rf <- randomForest(train.rdm, label.train.rf)
    
    pred.sam.rf <- predict(model.sam.rf, test.sam, type = "prob")[2]
    pred.mrmr.rf <- predict(model.mrmr.rf, test.mrmr, type = "prob")[2]
    pred.geode.rf <- predict(model.geode.rf, test.geode, type = "prob")[2]
    pred.rdm.rf <- predict(model.rdm.rf, test.rdm, type = "prob")[2]
    
    pred.sam.rf.list[i] <- pred.sam.rf
    pred.mrmr.rf.list[i] <- pred.mrmr.rf
    pred.geode.rf.list[i] <- pred.geode.rf
    pred.rdm.rf.list[i] <- pred.rdm.rf
    
    # use SVM to test the performance of three feature subsets
    train.sam.svm <- data.frame(cbind(train.sam, label.train))
    train.mrmr.svm <- data.frame(cbind(train.mrmr, label.train))
    train.geode.svm <- data.frame(cbind(train.geode, label.train))
    train.rdm.svm <- data.frame(cbind(train.rdm, label.train))
    
    test.sam.svm <- test.sam
    test.mrmr.svm <- test.mrmr
    test.geode.svm <- test.geode
    test.rdm.svm <- test.rdm
    
    model.sam.svm <- svm(formula = factor(label.train)~., data = train.sam.svm, probability = TRUE)
    model.mrmr.svm <- svm(formula = factor(label.train)~., data = train.mrmr.svm, probability = TRUE)
    model.geode.svm <- svm(formula = factor(label.train)~., data = train.geode.svm, probability = TRUE)
    model.rdm.svm <- svm(formula = factor(label.train)~., data = train.rdm.svm, probability = TRUE)
    
    pred.sam.svm <- predict(model.sam.svm, test.sam.svm, probability = TRUE)
    pred.mrmr.svm <- predict(model.mrmr.svm, test.mrmr.svm, probability = TRUE)
    pred.geode.svm <- predict(model.geode.svm, test.geode.svm, probability = TRUE)
    pred.rdm.svm <- predict(model.rdm.svm, test.rdm.svm, probability = TRUE)
    
    pred.sam.svm.list[i] <- attr(pred.sam.svm, "probabilities")[2]
    pred.mrmr.svm.list[i] <- attr(pred.mrmr.svm, "probabilities")[2]
    pred.geode.svm.list[i] <- attr(pred.geode.svm, "probabilities")[2]
    pred.rdm.svm.list[i] <- attr(pred.rdm.svm, "probabilities")[2]
    
    # RIDGE REGRESSION 
    set.seed(1234)
    ridge.cv.sam <- cv.glmnet(train.sam, label.train, alpha = 0, type.measure="mse",
                              family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    ridge.cv.mrmr <- cv.glmnet(train.mrmr, label.train, alpha = 0, type.measure="mse",
                               family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    ridge.cv.geode <- cv.glmnet(train.geode, label.train, alpha = 0, type.measure="mse",
                                family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    ridge.cv.rdm <- cv.glmnet(train.rdm, label.train, alpha = 0, type.measure="mse",
                              family = "binomial", standardize=F, nfolds = 3)
    
     # Optimal lambda:
    ridge.lam.sam <- ridge.cv.sam$lambda.min
    ridge.lam.mrmr <- ridge.cv.mrmr$lambda.min
    ridge.lam.geode <- ridge.cv.geode$lambda.min
    ridge.lam.rdm <- ridge.cv.rdm$lambda.min
    
    # Determine predicted values
    pred.sam.ridge <- predict(ridge.cv.sam, s = ridge.lam.sam, type="response", newx = test.sam)
    pred.mrmr.ridge <- predict(ridge.cv.mrmr, s = ridge.lam.mrmr, type="response", newx = test.mrmr)
    pred.geode.ridge <- predict(ridge.cv.geode, s = ridge.lam.geode, type="response", newx = test.geode)
    pred.rdm.ridge <- predict(ridge.cv.rdm, s = ridge.lam.rdm, type="response", newx = test.rdm)
    
    pred.sam.ridge.list[i] <- pred.sam.ridge
    pred.mrmr.ridge.list[i] <- pred.mrmr.ridge
    pred.geode.ridge.list[i] <- pred.geode.ridge
    pred.rdm.ridge.list[i] <- pred.rdm.ridge
    
    # LASSO REGRESSION 
    set.seed(1234)
    lasso.cv.sam <- cv.glmnet(train.sam, label.train, alpha = 1, type.measure="mse",
                              family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    lasso.cv.mrmr <- cv.glmnet(train.mrmr, label.train, alpha = 1, type.measure="mse",
                               family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    lasso.cv.geode <- cv.glmnet(train.geode, label.train, alpha = 1, type.measure="mse",
                                family = "binomial", standardize=F, nfolds = 3)
    set.seed(1234)
    lasso.cv.rdm <- cv.glmnet(train.rdm, label.train, alpha = 1, type.measure="mse",
                              family = "binomial", standardize=F, nfolds = 3)
    
     # Optimal lambda:
    lasso.lam.sam <- lasso.cv.sam$lambda.min
    lasso.lam.mrmr <- lasso.cv.mrmr$lambda.min
    lasso.lam.geode <- lasso.cv.geode$lambda.min
    lasso.lam.rdm <- lasso.cv.rdm$lambda.min
    
    # Determine predicted values
    pred.sam.lasso <- predict(lasso.cv.sam, s = lasso.lam.sam, type="response", newx = test.sam)
    pred.mrmr.lasso <- predict(lasso.cv.mrmr, s = lasso.lam.mrmr, type="response", newx = test.mrmr)
    pred.geode.lasso <- predict(lasso.cv.geode, s = lasso.lam.geode, type="response", newx = test.geode)
    pred.rdm.lasso <- predict(lasso.cv.rdm, s = lasso.lam.rdm, type="response", newx = test.rdm)
    
    pred.sam.lasso.list[i] <- pred.sam.lasso
    pred.mrmr.lasso.list[i] <- pred.mrmr.lasso
    pred.geode.lasso.list[i] <- pred.geode.lasso
    pred.rdm.lasso.list[i] <- pred.rdm.lasso
    
  }
  
  # calculate the stability
  stab.sam <- stab(select.genes.sam.list[[count.loop]])
  stab.mrmr <- stab(select.genes.mrmr.list[[count.loop]])
  stab.geode <- stab(select.genes.geode.list[[count.loop]])
  stab.rdm <- stab(select.genes.rdm.list[[count.loop]])
  
  stab.sam.list[count.loop] <- stab.sam
  stab.mrmr.list[count.loop] <- stab.mrmr
  stab.geode.list[count.loop] <- stab.geode
  stab.rdm.list[count.loop] <- stab.rdm
  
  # calculate the AUC
  auc.sam.rf <- auc(y.raw, pred.sam.rf.list)
  auc.mrmr.rf <- auc(y.raw, pred.mrmr.rf.list)
  auc.geode.rf <- auc(y.raw, pred.geode.rf.list)
  auc.rdm.rf <- auc(y.raw, pred.rdm.rf.list)

  auc.sam.svm <- auc(y.raw, pred.sam.svm.list)
  auc.mrmr.svm <- auc(y.raw, pred.mrmr.svm.list)
  auc.geode.svm <- auc(y.raw, pred.geode.svm.list)
  auc.rdm.svm <- auc(y.raw, pred.rdm.svm.list)
  
  auc.sam.ridge <- auc(y.raw, pred.sam.ridge.list)
  auc.mrmr.ridge <- auc(y.raw, pred.mrmr.ridge.list)
  auc.geode.ridge <- auc(y.raw, pred.geode.ridge.list)
  auc.rdm.ridge <- auc(y.raw, pred.rdm.ridge.list)
  
  auc.sam.lasso <- auc(y.raw, pred.sam.lasso.list)
  auc.mrmr.lasso <- auc(y.raw, pred.mrmr.lasso.list)
  auc.geode.lasso <- auc(y.raw, pred.geode.lasso.list)
  auc.rdm.lasso <- auc(y.raw, pred.rdm.lasso.list)
  
  auc.sam.rf.list[count.loop] <- auc.sam.rf
  auc.mrmr.rf.list[count.loop] <- auc.mrmr.rf
  auc.geode.rf.list[count.loop] <- auc.geode.rf
  auc.rdm.rf.list[count.loop] <- auc.rdm.rf
  
  auc.sam.svm.list[count.loop] <- auc.sam.svm
  auc.mrmr.svm.list[count.loop] <- auc.mrmr.svm
  auc.geode.svm.list[count.loop] <- auc.geode.svm
  auc.rdm.svm.list[count.loop] <- auc.rdm.svm
  
  auc.sam.ridge.list[count.loop] <- auc.sam.ridge
  auc.mrmr.ridge.list[count.loop] <- auc.mrmr.ridge
  auc.geode.ridge.list[count.loop] <- auc.geode.ridge
  auc.rdm.ridge.list[count.loop] <- auc.rdm.ridge
  
  auc.sam.lasso.list[count.loop] <- auc.sam.lasso
  auc.mrmr.lasso.list[count.loop] <- auc.mrmr.lasso
  auc.geode.lasso.list[count.loop] <- auc.geode.lasso
  auc.rdm.lasso.list[count.loop] <- auc.rdm.lasso
  
  print(paste('finished loop of', count.loop, '/', num.auc))
}

# plot the results
fig.auc.rf <- data.frame(num_select_genes = seq(num.select.genes.down, num.select.genes.up, step), sam.rf = auc.sam.rf.list, mrmr.rf = auc.mrmr.rf.list, geode.rf = auc.geode.rf.list, rdm.rf = auc.rdm.rf.list)
fig.auc.svm <- data.frame(num_select_genes = seq(num.select.genes.down, num.select.genes.up, step), sam.svm = auc.sam.svm.list, mrmr.svm = auc.mrmr.svm.list, geode.svm = auc.geode.svm.list, rdm.svm = auc.rdm.svm.list)
fig.auc.ridge <- data.frame(num_select_genes = seq(num.select.genes.down, num.select.genes.up, step), sam.ridge = auc.sam.ridge.list, mrmr.ridge = auc.mrmr.ridge.list, geode.ridge = auc.geode.ridge.list, rdm.ridge = auc.rdm.ridge.list)
fig.auc.lasso <- data.frame(num_select_genes = seq(num.select.genes.down, num.select.genes.up, step), sam.lasso = auc.sam.lasso.list, mrmr.lasso = auc.mrmr.lasso.list, geode.lasso = auc.geode.lasso.list, rdm.lasso = auc.rdm.lasso.list)
fig.auc.rf.long<-melt(fig.auc.rf,id='num_select_genes')
fig.auc.svm.long<-melt(fig.auc.svm,id='num_select_genes')
fig.auc.ridge.long<-melt(fig.auc.ridge,id='num_select_genes')
fig.auc.lasso.long<-melt(fig.auc.lasso,id='num_select_genes')
p1 <- ggplot(data=fig.auc.rf.long, aes(x=num_select_genes, y=value, colour=variable)) + geom_line() + xlab('Number of Selected Genes') + ylab('AUC of RF')
p2 <- ggplot(data=fig.auc.svm.long, aes(x=num_select_genes, y=value, colour=variable)) + geom_line() + xlab('Number of Selected Genes') + ylab('AUC of SVM')
p3 <- ggplot(data=fig.auc.ridge.long, aes(x=num_select_genes, y=value, colour=variable)) + geom_line() + xlab('Number of Selected Genes') + ylab('AUC of RIDGE')
p4 <- ggplot(data=fig.auc.lasso.long, aes(x=num_select_genes, y=value, colour=variable)) + geom_line() + xlab('Number of Selected Genes') + ylab('AUC of LASSO')
figure.auc <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
 
# # Fin
