#' Dataset containing an example LD profile
#'
#' A simulated LD profile, containing example LD statistics for
#' genetic distances of 0 to 0.0049, in bins of size 0.0001.
#'
#' @docType data
#'
#' @usage data(LDprofile)
#'
#' @format A data frame with 50 rows and 5 variables:
#' \describe{
#'   \item{bin}{the lower bound of each bin}
#'   \item{rsq}{the expected \eqn{r^2}{r^2} value for a pair of SNPs, where the genetic distance between them falls in the given bin}
#'   \item{sd}{the standard deviation of the expected \eqn{r^2}{r^2} value}
#'   \item{Beta_a}{the first shape parameter for the Beta distribution fitted for this bin}
#'   \item{Beta_b}{the second shape parameter for the Beta distribution fitted for this bin}
#' }
"LDprofile"
