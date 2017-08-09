\name{DensParcorr}
\alias{DensParcorr}
\title{
Conduct the Dens-Based approach for partial correlation estimation 
}
\description{
This function is to conduct the \emph{Dens}-based approach for partial correlation estimation in large scale brain network study.\cr
DensParcorr is the main function in this package. prec2dens and prec2part are sub-functions called by DensParcorr.
}
\usage{
DensParcorr(data,select=FALSE,dens.level="plateau",plateau.thresh=0.01,
            Parcorr.est=NULL,directory=NULL)
}

\arguments{
  \item{data}{Input data matrix with dimension of TxM where T is the number of observations and M is the number of nodes. For example, in fMRI data the T is the number of scans.}
  \item{select}{Whether to conduct the \emph{Dens}-based selection. If FALSE, output will only contain the estimated partial correlation list and precision matrix list corresponding to the default tuning parameter series ranging from 1e-8 to 0.6. If TRUE, the ouput will include the previous results and the selected partial correlation matrix and percision matrix corresponding to the specified density level. Default is FALSE. }
  \item{dens.level}{Specify the density level in \emph{Dens}-based tuning parameter selection method, including the plateau based density selection (\code{dens.level} = "plateau") and p percentage density selection (\code{dens.level} = p, 0<p<1). This option is valid only when \code{select}=TRUE. See Details. Default is "plateau".}
  \item{plateau.thresh}{The criterion to select the plateau. This option is valid only when \code{select}=TRUE. See Details. Default value is 0.01.}
  \item{Parcorr.est}{Previous output from DensParcorr function. }
  \item{directory}{The directory to output the figures and precision matrices and the partial correlation matrices. The default (\code{directory}=\code{NULL}) is to output in the current working directory.}
}

\details{
  This function implements the statistical method proposed in Wang et al. (2016) to estimate partial correlation matrix for studying direct connectivity in large-scale brain network. The method derives partial correlation based on the precision matrix estimated via Constrained L1-minimization Approach (CLIME) (Cai et al., 2011). This function applies the \emph{Dens}-based tuning parameter selection method in Wang et al. (2016) to help select an appropriate tuning parameter for sparsity control in the network estimation. Below is the breif step of \emph{Dens}-based approach.\cr
  
  First, we specify a series of tuning parameters \{\eqn{\lambda_n}\}. Then, based on \{\eqn{\lambda_n}\} we esimate a list of precision matices \eqn{\Omega(\lambda_n)} and and evaluate the density level of each precision matrix based on the \emph{Dens} criterion function in equation (5) of Wang et al. (2016). This will provide the users the profile of the density level corresponding to the series of tuning parameters in \{\eqn{\lambda_n}\}. Users can use the \code{dens.level} option to specify the desired density level in the precision matrix estimation. If \code{dens.level}="plateau", the function will select the plateau point \eqn{\lambda_{platu}} in the density profile based on the \code{plateau.thresh} and output the precision matrix \eqn{\Omega(\lambda_{platu})}. If \code{dens.level}=p and 0<p<1, the function will select the tuning parameter \eqn{\lambda_p} to achieve p percentage density and output the precision matrix \eqn{\Omega(\lambda_p)}. Then, the partial correlation matrix will be derived from the precision matrix. Further details can be found in the Reference. \cr
  
  The density profile and the heatmaps of precision matrices and partial correlation matrices will be saved in \code{directory}, and the esimated list of precision matrices and partial correlation matrices will also be saved in \code{directory}. \cr 
  
  When users would like to run the function multiple times on the same input data for different \code{dens.level}, it is computationally more efficient to read in the previous output from DensParcorr to \code{Parcorr.est} so that the function won't need to re-estimate the partial correlations based on the previous tuning parameters. 
}

\value{
  An R list from DensParcorr containing the following terms:
   \item{selected.partial.corr}{Selected Partial Correlation matrix corresponding to \code{dens.level}. Only when \code{select}=TRUE.}
   \item{selected.precision}{Selected Precision matrix corresponding to \code{dens.level}. Only when \code{select}=TRUE.}
   \item{selected.lambda}{Selected tuning parameter corresponding to \code{dens.level}. Only when \code{select}=TRUE.}
   \item{lambda.list}{The series of tuning parameters used for esimation and density profile.}
   \item{partial.corr.list}{Estimated Partial Correlation matrix list corresponding to \code{lambda.list}.}
   \item{precision.list}{Estimated Precision matrix list corresponding to \code{lambda.list}.}
   \item{Dens}{Actual density levels for estimated precision matrix list.}
   \item{Dens.Percentage}{Actual percentage density levels for estimated precision matrix list.}   
   \item{selection.method}{The method used for tuning parameter selection. For percentage \emph{Dens} selection, this value will include the actual \emph{Dens} precentage and the nominal \emph{Dens} percentage. Only when \code{select}=TRUE.}

  }

\references{
Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).   \emph{
     An efficient and reliable statistical method for estimating functional connectivity in large scale brain networks using partial correlation.
  } Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123 
  
Cai, T.T., Liu, W., and Luo, X. (2011).   \emph{
     A constrained \eqn{\ell_1}
  minimization approach for sparse precision matrix estimation.
  } Journal of the American Statistical Association 106(494): 594-607.
}

\author{Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
  Maintainer: Yikai Wang \email{yikai.wang@emory.edu} }
  
\examples{

# require(gplots)
# require(clime)

## Simulated the data to use. 
# data = matrix(rnorm(200),ncol=20)

##### Example 1: Estimate the partial correlation matrices for the 
##### default series of tuning parameters. 
# t0 = proc.time()[3]
# dens.est = DensParcorr(data,select=FALSE)
# proc.time()[3]-t0

##### Example 2: Estimate the network that reaches 40% density level. 
# partial.dens.est = DensParcorr(data,dens.level  =.4,select=TRUE)

###### Example 3: Now, estimate the 60% density level network based 
###### on the same data. To speed up computation, we read in the 
###### previous output from Example 2 into Parcorr.est
# t0 = proc.time()[3]
# partial.dens.est2 = DensParcorr(data, Parcorr.est = partial.dens.est,
#                                 dens.level=.6,select=TRUE)
# proc.time()[3]-t0

}