#' @title Conduct the Dens-Based approach for partial correlation estimation
#'
#' @description This function is to conduct the \emph{Dens}-based approach for partial
#' correlation estimation in large scale brain network study.\cr DensParcorr is
#' the main function in this package. prec2dens and prec2part are sub-functions
#' called by DensParcorr.
#'
#' @details This function implements the statistical method proposed in Wang et al.
#' (2016) to estimate partial correlation matrix for studying direct
#' connectivity in large-scale brain network. The method derives partial
#' correlation based on the precision matrix estimated via Constrained
#' L1-minimization Approach (CLIME) (Cai et al., 2011). This function applies
#' the \emph{Dens}-based tuning parameter selection method in Wang et al.
#' (2016) to help select an appropriate tuning parameter for sparsity control
#' in the network estimation. Below is the breif step of \emph{Dens}-based
#' approach.\cr
#'
#' First, we specify a series of tuning parameters {\eqn{\lambda_n}}. Then,
#' based on {\eqn{\lambda_n}} we esimate a list of precision matices
#' \eqn{\Omega(\lambda_n)} and and evaluate the density level of each precision
#' matrix based on the \emph{Dens} criterion function in equation (5) of Wang
#' et al. (2016). This will provide the users the profile of the density level
#' corresponding to the series of tuning parameters in {\eqn{\lambda_n}}. Users
#' can use the \code{dens.level} option to specify the desired density level in
#' the precision matrix estimation. If \code{dens.level}="plateau", the
#' function will select the plateau point \eqn{\lambda_{platu}} in the density
#' profile based on the \code{plateau.thresh} and output the precision matrix
#' \eqn{\Omega(\lambda_{platu})}. If \code{dens.level}=p and 0<p<1, the
#' function will select the tuning parameter \eqn{\lambda_p} to achieve p
#' percentage density and output the precision matrix \eqn{\Omega(\lambda_p)}.
#' Then, the partial correlation matrix will be derived from the precision
#' matrix. Further details can be found in the Reference. \cr
#'
#' The density profile and the heatmaps of precision matrices and partial
#' correlation matrices will be saved in \code{directory}, and the esimated
#' list of precision matrices and partial correlation matrices will also be
#' saved in \code{directory}. \cr
#'
#' When users would like to run the function multiple times on the same input
#' data for different \code{dens.level}, it is computationally more efficient
#' to read in the previous output from DensParcorr to \code{Parcorr.est} so
#' that the function won't need to re-estimate the partial correlations based
#' on the previous tuning parameters.
#'
#' @param data Input data matrix with dimension of TxM where T is the number of
#' observations and M is the number of nodes. For example, in fMRI data the T
#' is the number of scans.
#' @param select Whether to conduct the \emph{Dens}-based selection. If FALSE,
#' output will only contain the estimated partial correlation list and
#' precision matrix list corresponding to the default tuning parameter series
#' ranging from 1e-8 to 0.6. If TRUE, the ouput will include the previous
#' results and the selected partial correlation matrix and percision matrix
#' corresponding to the specified density level. Default is FALSE.
#' @param dens.level Specify the density level in \emph{Dens}-based tuning
#' parameter selection method, including the plateau based density selection
#' (\code{dens.level} = "plateau") and p percentage density selection
#' (\code{dens.level} = p, 0<p<1). This option is valid only when
#' \code{select}=TRUE. See Details. Default is "plateau".
#' @param plateau.thresh The criterion to select the plateau. This option is
#' valid only when \code{select}=TRUE. See Details. Default value is 0.01.
#' @param Parcorr.est Previous output from DensParcorr function.
#' @param directory The directory to output the figures and precision matrices
#' and the partial correlation matrices. The default
#' (\code{directory}=\code{NULL}) is to output in the current working
#' directory.
#'
#' @return An R list from DensParcorr containing the following terms:
#' \item{selected.partial.corr}{Selected Partial Correlation matrix
#' corresponding to \code{dens.level}. Only when \code{select}=TRUE.}
#' \item{selected.precision}{Selected Precision matrix corresponding to
#' \code{dens.level}. Only when \code{select}=TRUE.}
#' \item{selected.lambda}{Selected tuning parameter corresponding to
#' \code{dens.level}. Only when \code{select}=TRUE.} \item{lambda.list}{The
#' series of tuning parameters used for esimation and density profile.}
#' \item{partial.corr.list}{Estimated Partial Correlation matrix list
#' corresponding to \code{lambda.list}.} \item{precision.list}{Estimated
#' Precision matrix list corresponding to \code{lambda.list}.}
#' \item{Dens}{Actual density levels for estimated precision matrix list.}
#' \item{Dens.Percentage}{Actual percentage density levels for estimated
#' precision matrix list.} \item{selection.method}{The method used for tuning
#' parameter selection. For percentage \emph{Dens} selection, this value will
#' include the actual \emph{Dens} precentage and the nominal \emph{Dens}
#' percentage. Only when \code{select}=TRUE.}
#'
#' @author Yikai Wang, Jian Kang, Phebe Brenne Kemmer and Ying Guo\cr
#' Maintainer: Yikai Wang \email{yikai.wang@@emory.edu}
#' @references Wang, Y., Kang, J., Kemmer P. and Guo, Y. (2016).  \emph{ An
#' efficient and reliable statistical method for estimating functional
#' connectivity in large scale brain networks using partial correlation.  }
#' Front. Neurosci. 10:123. doi: 10.3389/fnins.2016.00123
#'
#' Cai, T.T., Liu, W., and Luo, X. (2011).  \emph{ A constrained \eqn{\ell_1}
#' minimization approach for sparse precision matrix estimation.  } Journal of
#' the American Statistical Association 106(494): 594-607.
#' @examples
#'
#'
#' # require(gplots)
#' # require(clime)
#'
#' ## Simulated the data to use.
#' data = matrix(rnorm(200),ncol=20)
#'
#' ##### Example 1: Estimate the partial correlation matrices for the
#' ##### default series of tuning parameters.
#' t0 = proc.time()[3]
#' dens.est = DensParcorr(data,select=FALSE, directory = tempfile())
#' proc.time()[3]-t0
#'
#' ##### Example 2: Estimate the network that reaches 40% density level.
#' partial.dens.est = DensParcorr(data,dens.level  =.4,select=TRUE,
#' directory = tempfile())
#'
#' ###### Example 3: Now, estimate the 60% density level network based
#' ###### on the same data. To speed up computation, we read in the
#' ###### previous output from Example 2 into Parcorr.est
#' t0 = proc.time()[3]
#' partial.dens.est2 = DensParcorr(data, Parcorr.est = partial.dens.est,
#'                                  dens.level=.6,select=TRUE,
#'                                  directory = tempfile())
#' proc.time()[3]-t0
#'
#' @importFrom clime clime
#' @import gplots
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend par plot points
#' @importFrom utils write.table
#' @export
DensParcorr <- function(
  data,
  select = FALSE,
  dens.level = "plateau",
  plateau.thresh = 0.01,
  Parcorr.est = NULL,
  directory = NULL
)
{
  if(is.null(directory)) {
    directory = file.path(getwd(), "DensParcorr.output")
  }

  lambda.max=0.6
  lambda.min=1e-8
  nlambda=10
  lambda = 10^(seq(log10(lambda.min), log10(lambda.max),length.out = nlambda))

  #### Calculate Precision Matrix ####
  data = as.matrix(data)

  if(!is.null(Parcorr.est))
  {
    if(is.null(Parcorr.est$precision.list)|is.null(Parcorr.est$lambda))
      stop("Parcorr.est is modified!")

    Prec.mat = Parcorr.est$precision.list
    lambda = Parcorr.est$lambda

    Prec.mat = Prec.mat[order(lambda)]
    lambda = lambda[order(lambda)]

  } else {
    lambda = lambda[order(lambda)]
    Prec.mat = clime::clime(data,lambda = lambda)$Omegalist
  }

  # Calculate Dens Values for Precision Matrix
  dens = vector()
  for(i in 1:length(Prec.mat))
  {
    dens[i]=prec2dens(Prec.mat[[i]])
  }

  # To guarantee the maximum of Dens is valid
  while(abs(1 - sort(dens)[length(dens)-1]/max(dens))>0.05)
  {
    lambda = c(min(lambda)/10,lambda)
    Prec.mat = append(clime(data,lambda = lambda[1])$Omegalist,Prec.mat)
    dens = c(prec2dens(Prec.mat[[1]]),dens)
  }

  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  #### Based on different Tuning Parameter Selection Method ####
  if(select)
  {

    if(dens.level=="plateau")
    {
      select.index = max(which((1-dens/max(dens))<=plateau.thresh))

      while(dens[select.index]==max(dens))
      {
        lam.min = lambda[select.index]
        lam.max = lambda[select.index+1]
        lambda2 = 10^(seq(log10(lam.min), log10(lam.max),length.out = 4))[2:3]
        Prec.mat2 = clime(data,lambda = lambda2)$Omegalist
        Prec.mat = append(Prec.mat2,Prec.mat)
        lambda = c(lambda2,lambda)

        Prec.mat = Prec.mat[order(lambda)]
        lambda = lambda[order(lambda)]
        for(i in 1:length(Prec.mat))
        {
          dens[i]=prec2dens(Prec.mat[[i]])
        }
        select.index = max(which((1-dens/max(dens))<=plateau.thresh))
      }

    }else if(is.numeric(dens.level)&dens.level<1&dens.level>0)
    {
      select.index = which.min(abs(dens-max(dens)*dens.level))

      while(abs(dens.level - dens[select.index]/max(dens))>0.05)
      {
        if(max(dens)*dens.level > dens[select.index])
        {
          lam.min = lambda[select.index-1]
          lam.max = lambda[select.index]
        }else
        {
          lam.min = lambda[select.index]
          lam.max = lambda[select.index+1]
        }
        lambda2 = 10^(seq(log10(lam.min), log10(lam.max),length.out = 4))[2:3]

        Prec.mat2 = clime(data,lambda = lambda2)$Omegalist
        Prec.mat = append(Prec.mat2,Prec.mat)
        lambda = c(lambda2,lambda)

        Prec.mat = Prec.mat[order(lambda)]
        lambda = lambda[order(lambda)]
        for(i in 1:length(Prec.mat))
        {
          dens[i]=prec2dens(Prec.mat[[i]])
        }
        select.index = which.min(abs(dens-max(dens)*dens.level))
      }
    }

    log10lam = -log(lambda)/log(10)



    png(filename = paste(directory,"/Dens.trace.plot.png",sep=""))
    par(mfrow=c(1,1))
    plot(log10lam[order(log10lam)],dens[order(log10lam)]/max(dens),xlab="-log10(lambda)",ylab="Percentage Dens",type="b")
    points(-log(lambda[select.index])/log(10),(dens/max(dens))[select.index],col="red",pch=19)
    legend("bottomright",legend="Selected",col="red",pch=19)
    dev.off()

  }else
  {
    png(filename = paste(directory,"/Dens.trace.plot.png",sep=""))
    par(mfrow=c(1,1))
    plot(-log(lambda[order(-log(lambda)/log(10))])/log(10),dens[order(-log(lambda)/log(10))]/max(dens),xlab="-log10(lambda)",ylab="Percentage Dens",type="b")
    dev.off()

  }

  #### Calculate Partial Correlation Matri from Precision Matrix ####

  dir.create(paste(directory,"/partial.correlation.matrix",sep=""))
  dir.create(paste(directory,"/precision.matrix",sep=""))

  Par.mat = list()
  for(i in 1:length(lambda))
  {
    Par.mat[[i]] = prec2part(Prec.mat[[i]])

    png(filename = paste(directory,"/partial.correlation.matrix/",i,".png",sep=""))
    heatmap.2(Par.mat[[i]],Rowv=F,Colv=F,scale="none",trace="none",density.info="none",xlab="",ylab="",
              main=paste("Dens Percentage=",round(dens[i]/max(dens),3),sep=""),col=bluered,dendrogram="none")
    dev.off()
    write.table(Par.mat[[i]],file =
                  paste(directory,"/partial.correlation.matrix/",i,"_Partial.DensPercentage_",
                        round(dens[i]/max(dens),3),".txt",sep=""),
                col.names = F,row.names = F)

    png(filename = paste(directory,"/precision.matrix/",i,".png",sep=""))
    heatmap.2(Par.mat[[i]],Rowv=F,Colv=F,scale="none",trace="none",density.info="none",xlab="",ylab="",
              main=paste("Dens Percentage=",round(dens[i]/max(dens),3),sep=""),col=bluered,dendrogram="none")
    dev.off()
    write.table(Prec.mat[[i]],file =
                  paste(directory,"/precision.matrix/",i,"_Precision.DensPercentage_",
                        round(dens[i]/max(dens),3),".txt",sep=""),
                col.names = F,row.names = F)
  }

  #### Summarize the results ####
  Results = list()
  if(select)
  {
    Results$selected.partial.corr = Par.mat[[select.index]]
    Results$selected.precision = Prec.mat[[select.index]]
    Results$selected.lambda = lambda[[select.index]]

    Results$partial.corr.list = Par.mat
    Results$precision.list = Prec.mat
    Results$lambda.list = lambda

  }else
  {
    Results$partial.corr.list = Par.mat
    Results$precision.list = Prec.mat
    Results$lambda.list = lambda
    Results$Dens = dens
    Results$Dens.Percentage = round(dens/max(dens),3)
  }

  if(select)
  {
    if(dens.level=="plateau")
    {
      Results$selection.method = "Dens-plateau"

    }else if(is.numeric(dens.level)&dens.level<1&dens.level>0)
    {
      Results$selection.method =  paste(round(dens.level*100,1),"% Dens (Actual=",round(dens[select.index]/max(dens)*100,1),"%)",sep="")
    }
    Results$Dens = dens
    Results$Dens.Percentage = round(dens/max(dens),3)
  }

  print("Figures, Estimated Precision Matrices and Partial Correlation are outputed in ")
  print(directory)

  return(Results)
}
