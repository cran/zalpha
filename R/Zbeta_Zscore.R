

#' Runs the Zbeta function using the Z score of the r-squared values for the region
#'
#' Returns a \eqn{Z_{\beta}^{Zscore}}{Zbeta} value for each SNP location supplied to the function, based on
#' the expected \eqn{r^2} values given an LD profile and genetic distances.
#' For more information about the \eqn{Z_{\beta}^{Zscore}}{Zbeta} statistic, please see Jacobs (2016).
#' The \eqn{Z_{\beta}^{Zscore}}{Zbeta} statistic is defined as:
#' \deqn{Z_{\beta}^{Zscore}=\frac{\sum_{i \in L,j \in R}\frac{r^2_{i,j}-E[r^2_{i,j}]}{\sigma[r^2_{i,j}]}}{|L||R|}}
#' where \code{|L|} and \code{|R|} are the number of SNPs to the left and right of the current locus within the given window \code{ws}, \eqn{r^2}{r^2} is equal to
#' the squared correlation between a pair of SNPs, \eqn{E[r^2]}{E[r^2]} is equal to the expected squared correlation between a pair of SNPs, given an LD profile, and \eqn{\sigma[r^2]}{\sigma[r^2]} is the standard deviation.
#'
#' The LD profile describes the expected correlation between SNPs at a given genetic distance, generated using simulations or
#' real data. Care should be taken to utilise an LD profile that is representative of the population in question. The LD
#' profile should consist of evenly sized bins of distances (for example 0.0001 cM per bin), where the value given is the (inclusive) lower
#' bound of the bin. Ideally, an LD profile would be generated using data from a null population with no selection, however one can be generated
#' using this data. See the \code{\link{create_LDprofile}} function for more information on how to create an LD profile.
#'
#' @importFrom stats cor na.omit
#'
#' @param pos A numeric vector of SNP locations
#' @param ws The window size which the \eqn{Z_{\beta}^{Zscore}}{Zbeta} statistic will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param x A matrix of SNP values. Columns represent chromosomes; rows are SNP locations. Hence, the number of rows should equal the length of the \code{pos} vector. SNPs should all be biallelic.
#' @param dist A numeric vector of genetic distances (e.g. cM, LDU). This should be the same length as \code{pos}.
#' @param LDprofile_bins A numeric vector containing the lower bound of the bins used in the LD profile. These should be of equal size.
#' @param LDprofile_rsq A numeric vector containing the expected \eqn{r^2}{r^2} values for the corresponding bin in the LD profile. Must be between 0 and 1.
#' @param LDprofile_sd A numeric vector containing the standard deviation of the \eqn{r^2}{r^2} values for the corresponding bin in the LD profile.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistic to be calculated. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate \eqn{Z_{\beta}^{Zscore}}{Zbeta} for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate \eqn{Z_{\beta}^{Zscore}}{Zbeta} for every SNP in the \code{pos} vector.
#'
#' @return A list containing the SNP positions and the \eqn{Z_{\beta}^{Zscore}}{Zbeta} values for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @examples
#' ## load the snps and LDprofile example datasets
#' data(snps)
#' data(LDprofile)
#' ## run Zbeta_Zscore over all the SNPs with a window size of 3000 bp
#' Zbeta_Zscore(snps$bp_positions,3000,as.matrix(snps[,3:12]),snps$cM_distances,
#'  LDprofile$bin,LDprofile$rsq,LDprofile$sd)
#' ## only return results for SNPs between locations 600 and 1500 bp
#' Zbeta_Zscore(snps$bp_positions,3000,as.matrix(snps[,3:12]),snps$cM_distances,
#'  LDprofile$bin,LDprofile$rsq,LDprofile$sd,X=c(600,1500))
#'
#' @export
#' @seealso \code{\link{create_LDprofile}}
Zbeta_Zscore<-function(pos, ws, x, dist, LDprofile_bins, LDprofile_rsq, LDprofile_sd, minRandL = 4, minRL = 25, X = NULL){

  #Check things are in the correct format

  #Check pos is a numeric vector
  if (is.numeric(pos) ==FALSE || is.vector(pos)==FALSE){
    stop("pos must be a numeric vector")
  }
  #Check x is a matrix
  if (is.matrix(x)==FALSE){
    stop("x must be a matrix")
  }
  #Check x has rows equal to the length of pos
  if (length(pos) != nrow(x)){
    stop("The number of rows in x must equal the number of SNP locations given in pos")
  }
  #Check SNPs are all biallelic
  if(sum(apply(x,1,function(x){length(na.omit(unique(x)))}) != 2)>0){
    stop("SNPs must all be biallelic")
  }
  #Check dist is a numeric vector
  if (is.numeric(dist) ==FALSE || is.vector(dist)==FALSE){
    stop("dist must be a numeric vector")
  }
  #Check dist is the same length as pos
  if (length(pos) != length(dist)){
    stop("The number of values in dist must equal the number of SNP locations given in pos")
  }
  #Check windowsize is a number greater than 0
  if(is.numeric(ws) ==FALSE || ws <= 0){
    stop("ws must be a number greater than 0")
  }
  #Check LDprofile_bins is a numeric vector
  if (is.numeric(LDprofile_bins) ==FALSE || is.vector(LDprofile_bins)==FALSE){
    stop("LDprofile_bins must be a numeric vector")
  }
  #Get bin size from LDprofile_bins
  bin_size<-LDprofile_bins[2]-LDprofile_bins[1]
  #Check LDprofile_bins are of equal size
  if (isTRUE(all.equal(diff(LDprofile_bins),rep(bin_size,length(LDprofile_bins)-1)))==FALSE){
    stop("LDprofile_bins must be of equal size")
  }
  #Check LDprofile_rsq is a numeric vector
  if (is.numeric(LDprofile_rsq) ==FALSE || is.vector(LDprofile_rsq)==FALSE){
    stop("LDprofile_rsq must be a numeric vector")
  }
  #Check LDprofile_sd is a numeric vector
  if (is.numeric(LDprofile_sd) ==FALSE || is.vector(LDprofile_sd)==FALSE){
    stop("LDprofile_sd must be a numeric vector")
  }
  #Check values of LDprofile_rsq are between 0 and 1
  if (sum(LDprofile_rsq<0 | LDprofile_rsq>1)>0){
    stop("Values stored in LDprofile_rsq must be between 0 and 1")
  }
  #Check that the LDprofile vectors are the same length
  if (length(LDprofile_bins) != length(LDprofile_rsq)){
    stop("LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_bins")
  }
  #Check that the LDprofile vectors are the same length
  if (length(LDprofile_bins) != length(LDprofile_sd)){
    stop("LDprofile_sd must contain the same number of values as there are bins given in LDprofile_bins")
  }
  #Check minRandL is 0 or greater
  if(is.numeric(minRandL) ==FALSE || minRandL < 0){
    stop("minRandL must be a number greater than or equal to 0")
  }
  #Check minRL is 0 or greater
  if(is.numeric(minRL) ==FALSE || minRL < 0){
    stop("minRL must be a number greater than or equal to 0")
  }
  #If X is specified, check it is in the correct format
  if (is.null(X)==FALSE){
    if(is.numeric(X)==FALSE || is.vector(X)==FALSE){
      stop("X should be a numeric vector of length 2 e.g. c(100,200)")
    } else {
      if (length(X) != 2){
        stop("X should be a numeric vector of length 2 e.g. c(100,200)")
      } else {
        # X is in the correct format
        # Check that X will actually return a result (i.e. that the region specied by X overlaps with pos)
        if ((length(pos[pos>=X[1] & pos <= X[2]])>0) == FALSE){
          stop("The region specified by X is outside the region contained in the pos vector")
        }
      }
    }
  } else {
    # Set X equal to the extremes of pos
    X<-c(pos[1],pos[length(pos)])
  }

  # Force the R code to print decimals in full rather than in scientific format
  oldOptions<-options(scipen=999)
  on.exit(options(oldOptions))

  #Change matrix x to numeric if it isn't already
  if (is.numeric(x)==FALSE){
    x<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
  }

  # Set up output list
  outputLength<-length(pos[pos>=X[1] & pos <= X[2]])
  outputList<-list(position=pos[pos>=X[1] & pos <= X[2]],Zbeta_Zscore=rep(NA,outputLength))


  # Loop over each position in the output list and calculate the expected Zbeta
  for (i in 1:outputLength){

    # Current physical position in chromosome
    currentPos<-outputList$position[i]

    ## check L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    if  (noL < minRandL || noR < minRandL || noL*noR < minRL){
      #NA
      outputList$Zbeta_Zscore[i]<-NA
    } else {
      # Find distances between each SNP in L and round to bin size
      bins<-sapply(outer(dist[pos<=currentPos+ws/2 & pos > currentPos],dist[pos>=currentPos-ws/2 & pos < currentPos],"-"),assign_bins,bin_size=bin_size)
      bins[bins>max(LDprofile_bins)]<-max(LDprofile_bins)
      rsq<-as.vector(t((cor(t(x[pos>=currentPos-ws/2 & pos<=currentPos+ws/2,]),use="pairwise.complete.obs")^2)[1:noL,(noL+2):(noL+noR+1)]))
      rsqExp<-merge(data.frame(bins=as.character(bins),rsq),data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_rsq,LDprofile_sd),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
      rsqSum<-sum((rsqExp$rsq-rsqExp$LDprofile_rsq)/rsqExp$LDprofile_sd)

      outputList$Zbeta_Zscore[i]<-rsqSum/(noL*noR)
    }
  }
  if (sum(is.na(outputList$Zbeta_Zscore))==outputLength){
    warning("No Zbeta_Zscore values were calculated, try reducing minRandL and minRL or increasing the window size")
  }
  return(outputList)
}

