#' Calculate Partial Correlation Matrix from Precision Matrix
#'
#' This function is to derive the partial correlation matrix from the precision
#' matrix.
#'
#' @param Precision Input precision matrix (symmetric and positive definite).
#' @return Calculated partial correlation matrix.
#' @author Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
#' Maintainer: Yikai Wang \email{yikai.wang@@emory.edu}
#' @references Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).  \emph{ An
#' efficient and reliable statistical method for estimating functional
#' connectivity in large scale brain networks using partial correlation.  }
#' Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123
#' @export
prec2part <- function(Precision)
{
  Precision = as.matrix(Precision)
  D.Prec = diag(diag(Precision)^(-.5))
  return(diag(2,dim(Precision)[1])-D.Prec%*%Precision%*%D.Prec)
}
