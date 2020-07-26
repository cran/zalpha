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
results<-Zalpha(snps$bp_positions,3000,as.matrix(snps[,3:12]))
results
plot(results$position,results$Zalpha)

## -----------------------------------------------------------------------------
Zalpha(snps$bp_positions,3000,as.matrix(snps[,3:12]),X=c(500,1000))

## -----------------------------------------------------------------------------
snps$cM_distances

## -----------------------------------------------------------------------------
data(LDprofile)
LDprofile

## -----------------------------------------------------------------------------
Zalpha_expected(snps$bp_positions, 3000, snps$cM_distances, LDprofile$bin, LDprofile$rsq)

## -----------------------------------------------------------------------------
Zalpha_rsq_over_expected(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, LDprofile$bin, LDprofile$rsq)
Zalpha_log_rsq_over_expected(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, LDprofile$bin, LDprofile$rsq)
Zalpha_Zscore(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, LDprofile$bin, LDprofile$rsq, LDprofile$sd)
Zalpha_BetaCDF(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, LDprofile$bin, LDprofile$Beta_a, LDprofile$Beta_b)

## -----------------------------------------------------------------------------
results<-Zbeta(snps$bp_positions,3000,as.matrix(snps[,3:12]))
results
plot(results$position,results$Zbeta)

## -----------------------------------------------------------------------------
Zbeta_expected(snps$bp_positions, 3000, snps$cM_distances,
               LDprofile$bin, LDprofile$rsq)
Zbeta_rsq_over_expected(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, 
                        LDprofile$bin, LDprofile$rsq)
Zbeta_log_rsq_over_expected(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, 
                            LDprofile$bin, LDprofile$rsq)
Zbeta_Zscore(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, 
             LDprofile$bin, LDprofile$rsq, LDprofile$sd)
Zbeta_BetaCDF(snps$bp_positions, 3000, as.matrix(snps[,3:12]), snps$cM_distances, 
              LDprofile$bin, LDprofile$Beta_a, LDprofile$Beta_b)

## -----------------------------------------------------------------------------
Zalpha_all(snps$bp_positions,3000,as.matrix(snps[,3:12]))

## -----------------------------------------------------------------------------
create_LDprofile(snps$cM_distances,as.matrix(snps[,3:12]),bin_size = 0.001,beta_params = TRUE)

## -----------------------------------------------------------------------------
## Generate three chromosomes of data - cM distances and SNP values
chrom1_cM_distances<-snps$cM_distances
chrom1_snp_values<-as.matrix(snps[,3:12])

chrom2_cM_distances<-snps$cM_distances
chrom2_snp_values<-as.matrix(snps[,3:12])

chrom3_cM_distances<-snps$cM_distances
chrom3_snp_values<-as.matrix(snps[,3:12])

## create a list of the cM distances
cM_distances_list<-list(chrom1_cM_distances,chrom2_cM_distances,chrom3_cM_distances)

## create a list of SNP value matrices
snp_values_list<-list(chrom1_snp_values,chrom2_snp_values,chrom3_snp_values)

## create the LD profile using the lists as the dist and x parameters
create_LDprofile(cM_distances_list,snp_values_list,bin_size = 0.001,beta_params = TRUE)

