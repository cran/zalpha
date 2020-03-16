## ---- echo = FALSE, message=FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(zalpha)

## -----------------------------------------------------------------------------
library(zalpha)
data(snps)
## This is what the dataset looks like:
snps

## -----------------------------------------------------------------------------
results<-Zalpha(snps$positions,3000,as.matrix(snps[,3:12]))
results
plot(results$position,results$Zalpha)

## -----------------------------------------------------------------------------
Zalpha(snps$positions,3000,as.matrix(snps[,3:12]),X=c(500,1000))

## -----------------------------------------------------------------------------
LDprofile

## -----------------------------------------------------------------------------
Zalpha_expected(snps$positions, 3000, snps$distances, LDprofile$bin, LDprofile$rsq)

## -----------------------------------------------------------------------------
Zalpha_all(snps$positions,3000,as.matrix(snps[,3:12]))

