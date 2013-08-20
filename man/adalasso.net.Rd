\name{adalasso.net}
\alias{adalasso.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Partial Correlations with (Adaptive) Lasso}
\description{
 This function computes the matrix of partial correlations based on an
  estimation of the corresponding regression models via lasso and adaptive lasso respectively.

}
\usage{
adalasso.net(X, k = 10,use.Gram=FALSE,both=TRUE,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables.}
  \item{k}{the number of splits in  \code{k}-fold cross-validation. The
    same \code{k} is used for the estimation of the weights and the
    estimation of the penalty term for adaptive lasso. Default value is \code{k}=10. }
    \item{use.Gram}{When the number of variables is very large, you may not want LARS to precompute the Gram matrix. Default is \code{use.Gram}=FALSE. }
    \item{both}{Logical. If both=FALSE, only the lasso solution is computed. Default is both=TRUE.}
    \item{verbose}{Print information on conflicting signs etc. Default is \code{verbose=FALSE}}
}
\details{
 For each of the columns of \code{X}, a regression model based on
  (adaptive) lasso is computed. In each of the \code{k}-fold cross-validation steps, the weights for adaptive lasso are computed in 
  terms of a lasso fit. (The optimal value of the
  penalty term is selected via \code{k}-fold cross-validation). Note that this implies that a lasso solution is computed k*k times! Finally, the results of the regression models are
  transformed via the function \code{Beta2parcor}.
}
\value{
  \item{pcor.adalasso}{estimated matrix of partial correlation coefficients for adaptive lasso.}
  \item{pcor.lasso}{estimated matrix of partial correlation coefficients for lasso.}
  ...
}
\references{H. Zou (2006) "The Adaptive Lasso and its Oracle Property",
  Journal of the American Statistical Association. 101 (476):  1418-1429.

  
  N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks using Gaussian Graphical Models", BMC Bioinformatics, 10:384

\url{http://www.biomedcentral.com/1471-2105/10/384/}
}

\author{Nicole Kraemer}

\seealso{\code{\link{Beta2parcor}}, \code{\link{adalasso}}}
\examples{
n<-20
p<-10
X<-matrix(rnorm(n*p),ncol=p)
pc<-adalasso.net(X,k=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
