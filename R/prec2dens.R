#' Calculate Dens value from Precision Matrix
#'
#' This function evaluates the density level of a precision matrix based on the
#' Dens criterion function in equation (5) of Wang et al. (2016).
#'
#'
#' @param precision Input precision matrix (symmetric and positive definite).
#' @return Density level from Precision matrix.
#' @author Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
#' Maintainer: Yikai Wang \email{yikai.wang@@emory.edu}
#' @references Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).  \emph{ An
#' efficient and reliable statistical method for estimating functional
#' connectivity in large scale brain networks using partial correlation.  }
#' Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123
#' @export
prec2dens<- function(precision)
{
  return(sum(abs(precision)))
}
