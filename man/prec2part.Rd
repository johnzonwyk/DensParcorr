\name{prec2part}
\alias{prec2part}
\title{
Calculate Partial Correlation Matrix from Precision Matrix
}
\description{
This function is to derive the partial correlation matrix from the precision matrix.
}
\usage{
prec2part(Precision)
}
\arguments{
  \item{Precision}{Input precision matrix (symmetric and positive definite).}
}
\value{
  Calculated partial correlation matrix.}

\references{Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).   \emph{
     An efficient and reliable statistical method for estimating functional connectivity in large scale brain networks using partial correlation.
  } Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123
}

\author{Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
  Maintainer: Yikai Wang \email{yikai.wang@emory.edu} }
