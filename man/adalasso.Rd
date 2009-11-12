\name{adalasso}
\alias{adalasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adaptive Lasso}
\description{
  This function computes the lasso and adaptive lasso solution based on k-fold cross-validation. The initial
  weights for adaptive lasso are computed from a lasso fit.

}
\usage{
adalasso(X, y, k = 10, use.Gram = TRUE,both=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables}
  \item{y}{vector of responses. The length of y must equal the number of
  rows of X}
  \item{k}{the number of splits in  \code{k}-fold cross-validation. The
    same \code{k} is used for the estimation of the weights and the
    estimation of the penalty term for adaptive lasso. Default is k=10.}

  \item{use.Gram}{When the number of variables is very large, you may not want LARS to precompute the Gram matrix. Default is \code{use.Gram}=TRUE. }
    \item{both}{Logical. If both=FALSE, only the lasso coefficients are computed. Default is both=TRUE.} 
}

\details{
  In each of the \code{k}-fold cross-validation steps, the weights for adaptive lasso are computed in 
  terms of a lasso fit. (The optimal value of the
  penalty term is selected via \code{k}-fold cross-validation). Note that this implies that a lasso solution is computed k*k times!
}
\value{
  \item{coefficients.adalasso}{regression coefficients for adaptive lasso.}
  \item{coefficients.lasso}{regression coefficients for lasso.}
  \item{cv.lasso}{cv error for the optimal lasso model.}
  \item{cv.adalasso}{cv error for the optimal adaptive lasso model.}
}
\references{H. Zou (2006) "The Adaptive Lasso and its Oracle Property",
  Journal of the American Statistical Association 101 (476):  1418-1429.

  
  N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks using Gaussian Graphical Models", BMC Bioinformatics, in press}

\author{Nicole Kraemer, Juliane Schaefer}

\seealso{\code{\link{Beta2parcor}}, \code{\link{adalasso.net}}}

\examples{
n<-100 # number of observations
p<-60 # number of variables
X<-matrix(rnorm(n*p),ncol=p) 
y<-rnorm(n)
ada.object<-adalasso(X,y,k=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
