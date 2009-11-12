\name{mylars}
\alias{mylars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for Lasso}
\description{
 This function computes the cross-validation-optimal regression coefficients for lasso.
}
\usage{
mylars(X, y, k = 10,use.Gram=TRUE,normalize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{matrix of observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables}
   \item{y}{vector of responses. The length of y must equal the number of
  rows of X}
  \item{k}{the number of splits in  \code{k}-fold cross-validation.  Default is k=10.}
    \item{use.Gram}{When the number of variables is very large, you may not want LARS to precompute the Gram matrix. Default is \code{use.Gram}=TRUE. }
    \item{normalize}{Should the columns of \code{X} be scaled? Default is normalize=TRUE.}
}


\details{
 We use the glmnet() function from the glmnet package to compute the fit. Note that in Kraemer et. al. (2009), we used the lars() 
 function from the lars package, which is much slower than glmnet().}


\value{
\item{coefficients}{cross-validation optimal regression coefficients, without intercept.}
\item{cv.lasso}{cv error for the optimal model.}
}
\references{R. Tibshirani (1997) "Regression Shrinkage and Selection via the
  Lasso", Journal of the Royal Statistical Society B, 58 (1)


 

N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", BMC Bioinformatics, in press

}
\author{Nicole Kraemer}

\seealso{\code{\link{Beta2parcor}}, \code{\link{adalasso}}}
\examples{

n<-20
p<-50
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)
dummy<-mylars(X,y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
