
#' Runs the Zalpha function
#'
#' Returns a \eqn{Z_{\alpha}}{Zalpha} value for each SNP location supplied to the function.
#' For more information about the \eqn{Z_{\alpha}}{Zalpha} statistic, please see Jacobs (2016).
#' The \eqn{Z_{\alpha}}{Zalpha} statistic is defined as:
#' \deqn{Z_{\alpha}=\frac{{|L| \choose 2}^{-1}\sum_{i,j \in L}r^2_{i,j} + {|R| \choose 2}^{-1}\sum_{i,j \in L}r^2_{i,j}}{2}}
#' where \code{|L|} and \code{|R|} are the number of SNPs to the left and right of the current locus within the given window \code{ws}, and \eqn{r^2}{r^2} is equal to the squared correlation between a pair of SNPs
#'
#' @importFrom stats cor na.omit
#'
#' @param pos A numeric vector of SNP locations
#' @param ws The window size which the \eqn{Z_{\alpha}}{Zalpha} statistic will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param x A matrix of SNP values. Columns represent chromosomes; rows are SNP locations. Hence, the number of rows should equal the length of the \code{pos} vector. SNPs should all be biallelic.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistic to be calculated. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate \eqn{Z_{\alpha}}{Zalpha} for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate \eqn{Z_{\alpha}}{Zalpha} for every SNP in the \code{pos} vector.
#'
#' @return A list containing the SNP positions and the \eqn{Z_{\alpha}}{Zalpha} values for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @examples
#' ## load the snps example dataset
#' data(snps)
#' ## run Zalpha over all the SNPs with a window size of 3000 bp
#' Zalpha(snps$bp_positions,3000,as.matrix(snps[,3:12]))
#' ## only return results for SNPs between locations 600 and 1500 bp
#' Zalpha(snps$bp_positions,3000,as.matrix(snps[,3:12]),X=c(600,1500))
#'
#' @export

Zalpha <- function(pos, ws, x, minRandL = 4, minRL = 25, X = NULL) {
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
  #Check windowsize is a number greater than 0
  if(is.numeric(ws) ==FALSE || ws <= 0){
    stop("ws must be a number greater than 0")
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

  #Change matrix x to numeric if it isn't already
  if (is.numeric(x)==FALSE){
    x<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
  }

  # Set up output list
  outputLength<-length(pos[pos>=X[1] & pos <= X[2]])
  outputList<-list(position=pos[pos>=X[1] & pos <= X[2]],Zalpha=rep(NA,outputLength))


  # Loop over each position in the output list and calculate Zalpha
  for (i in 1:outputLength){

    # Current physical position in chromosome
    currentPos<-outputList$position[i]

    ## check L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    if  (noL < minRandL || noR < minRandL || noL*noR < minRL){
      #NA
      outputList$Zalpha[i]<-NA
    } else {
      ##Left
      LrsqSum<-sum(lower_triangle(cor(t(x[pos>=currentPos-ws/2 & pos < currentPos,]),use="pairwise.complete.obs")^2))
      ##Right
      RrsqSum<-sum(lower_triangle(cor(t(x[pos<=currentPos+ws/2 & pos > currentPos,]),use="pairwise.complete.obs")^2))
      outputList$Zalpha[i]<-(LrsqSum/choose(noL,2)+RrsqSum/choose(noR,2))/2
    }
  }
  if (sum(is.na(outputList$Zalpha))==outputLength){
    warning("No Zalpha values were calculated, try reducing minRandL and minRL or increasing the window size")
  }
  return(outputList)
}
