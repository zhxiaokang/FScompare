library(openxlsx)
# prepare for the data set
data.raw <- read.xlsx("/Users/xzh004/Desktop/research/dCod/RNA-Seq_analysis/UiB_Goksoyr_dCod_1_2017/PCLS_BaP-EE2_RNAseq/oldTrans/countFile/normCount/TPM_Normalised_Counts_ALL.xlsx", sheetIndex = 1)  # each row is a protein, each column is a sample, but the first two columns are pretein names and id

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

# build the data frame
