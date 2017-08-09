\name{prec2dens}
\alias{prec2dens}
\title{
Calculate Dens value from Precision Matrix
}
\description{
This function evaluates the density level of a precision matrix based on the Dens criterion function in equation (5) of Wang et al. (2016).  
}
\usage{
prec2dens(precision)
}
\arguments{
  \item{precision}{Input precision matrix (symmetric and positive definite).}
}
\value{
  Density level from Precision matrix.}

\references{Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).   \emph{
     An efficient and reliable statistical method for estimating functional connectivity in large scale brain networks using partial correlation.
  } Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123
}

\author{Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
  Maintainer: Yikai Wang \email{yikai.wang@emory.edu} }
