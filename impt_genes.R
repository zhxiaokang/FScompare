# Select the features with different feature selection methods with all samples
# The purpose is to provide the "important" protein lists for downstream analysis

#remove(list = ls())

library(xlsx)
library(samr)
library(mRMRe)
library(pROC)

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

# prepare for the data set
data.raw <- read.xlsx("./data/PCB153.xlsx", "data_raw")  # each row is a protein, each column is a sample, but the first two columns are pretein names and id

# set the number of selected proteins that you want
num.select.genes <- 100

nsample <- 20  # number of samples
nfeature <- nrow(data.raw) - 1  # number of features (proteins)

x.raw <- log2(data.raw[ , c(3:12, 33:42)] + 1)   # nsample samples (control group and high dose group)
num.control <- 10
num.treat <- 10
if (nsample != num.control + num.treat){
  print('Number of samples do not agree!')
}
y.raw <- c(rep(1,num.control), rep(2,num.treat))  # control samples are labeled as 0, and treatment samples are labeled as 1
pr.name <- data.raw$accession
pr.id <- data.raw$id

# use SAM to select the significantly differentially regulated genes (proteins)
select.genes.sam.names <- rep(0, num.select.genes)
data <- list(x=x.raw,y=y.raw,genenames=paste("gene",as.character(pr.id),sep=""), geneid=as.character(pr.id), logged2=TRUE)  
invisible(capture.output(samr.obj <- samr(data, resp.type="Two class unpaired", nperms=100)))
invisible(capture.output(delta.table <- samr.compute.delta.table(samr.obj)))
del <- -3
siggenes.table <- samr.compute.siggenes.table(samr.obj, del, data, delta.table)
siggenes.table$genes.lo[,4] <- -as.numeric((siggenes.table$genes.lo[,4]))  # reverse the Score(d) of genes.lo
siggenes <- rbind(siggenes.table$genes.up, siggenes.table$genes.lo)  # combine the up and low genes
siggenes <- siggenes[order(siggenes[,4],decreasing = TRUE),]  # sort the genes aucording to the Score(d)
select.genes.sam <- siggenes[1:num.select.genes, 1]  # significantly differentially regulated genes selected by SAM
select.genes.sam <- as.numeric(select.genes.sam) - 1
for (s in seq(1, num.select.genes)){
  select.genes.sam.names[s] <- as.character(data.raw$accession[select.genes.sam[s]])
}

# use Feature Selection method mRMR to select the important genes (proteins)
select.genes.mrmr.names <- rep('protein_accession', num.select.genes)
data.mrmr <- cbind(y.raw, t(x.raw))  # in mRMRe, each row is a sample, and each colunm is a feature, and the first colunm is class label
dd <- mRMR.data(data = data.frame(data.mrmr))
result.mrmr <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 1, feature_count = num.select.genes)
select.genes.mrmr <- result.mrmr@filters
select.genes.mrmr <- as.numeric(select.genes.mrmr$`1`) - 1
for (s in seq(1, num.select.genes)){
  select.genes.mrmr.names[s] <- as.character(data.raw$accession[select.genes.mrmr[s]])
}

# use GeoDE to select the important genes (proteins)
select.genes.geode <- rep('protein_accession', num.select.genes)
gammas <- 1
data.geode <- cbind(pr.name, x.raw)
chdir_analysis <- chdirAnalysis(data.geode, factor(y.raw), gammas, CalculateSig=FALSE, nnull=10)
select.genes.geode.names <- names(chdir_analysis$results[[1]][1:num.select.genes])
for (s in seq(1, num.select.genes)){
  select.genes.geode[s] <- which(data.raw$accession == select.genes.geode.names[s])
}

# save the protein lists to txt files
write.table(select.genes.sam.names,'./output/pcb_high/protein_sam.txt',quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(select.genes.mrmr.names,'./output/pcb_high/protein_mrmr.txt',quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(select.genes.geode.names,'./output/pcb/high/protein_geode.txt',quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


