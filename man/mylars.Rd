\name{mylars}
\alias{mylars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for Lasso}
\description{
 This function computes the cross-validation-optimal regression coefficients for lasso.
}
\usage{
mylars(X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE)
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
     \item{intercept}{Should an intercept be included? Default is \code{intercept=TRUE}.}
}


\details{
 We use the glmnet() function from the glmnet package to compute the fit. Note that in Kraemer et. al. (2009), we used the lars() 
 function from the lars package, which is much slower than glmnet().}


\value{
\item{lambda}{vector of paramter values from which the optimal parameter is selected}
\item{cv}{cross-validated error for all \code{lambda} values}
\item{lambda.opt}{cross-validation optimal parameter}
\item{cv.lasso}{cv error for the optimal model.}
\item{intercept}{cross-validation optimal intercept. If \code{intercept=FALSE} was specified, the intercept is set to 0.}
\item{coefficients}{cross-validation optimal regression coefficients, without intercept}
}
\references{R. Tibshirani (1997) "Regression Shrinkage and Selection via the
  Lasso", Journal of the Royal Statistical Society B, 58 (1)


 

N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", BMC Bioinformatics, 10:384

\url{http://www.biomedcentral.com/1471-2105/10/384/}

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
