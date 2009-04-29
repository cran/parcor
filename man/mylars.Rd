\name{mylars}
\alias{mylars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for Lasso}
\description{
 This function computes the cross-validation-optimal regression coefficients for lasso.
}
\usage{
mylars(X, y, k = 10,fraction = seq(from = 0, to = 1, length = 1000),use.Gram=TRUE,normalize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{matrix of observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables}
   \item{y}{vector of responses. The length of y must equal the number of
  rows of X}
  \item{k}{the number of splits in  \code{k}-fold cross-validation.  Default is k=10.}
    \item{use.Gram}{When the number of variables is very large, you may not want LARS to precompute the Gram matrix. Default is \code{use.Gram}=TRUE. }
    \item{fraction}{vector of possible regularization parameters, in the range from 0 to 1.}
    \item{normalize}{Should the columns of \code{X} be scaled? Default is normalize=TRUE.}
}
\details{
 This is a variation of the \code{cv.lars} function of the \code{lars} package. Here, we adjust the regularization parameter \code{fraction} in order to avoid 
 its peaking behavior in the n=p case. See Kraemer (2009) for more details.}
\value{
\item{coefficients}{cross-validation optimal regression coefficients, without intercept.}
\item{cv.lasso}{cv error for the optimal model.}
}
\references{R. Tibshirani (1997) "Regression Shrinkage and Selection via the
  Lasso", Journal of the Royal Statistical Society B, 58 (1)

N. Kraemer (2009) "On the Peaking Phenomenon of the Lasso in Model Selection", preprint, 
\url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}
 

N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", preprint

  \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}

}
\author{Nicole Kraemer}

\seealso{\code{\link{Beta2parcor}}, \code{\link{adalasso}}}
\examples{

n<-50
p<-10
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)
dummy<-mylars(X,y,k=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
