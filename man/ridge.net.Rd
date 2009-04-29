\name{ridge.net}
\alias{ridge.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial correlations with ridge regression.}
\description{
  This function computes the matrix of partial correlations via an
  estimation of the corresponding regression models via Ridge Regression.
}
\usage{
ridge.net(X, lambda, plot.it = FALSE, scale = TRUE, k = 10,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables.}
  \item{lambda}{Vector of penalty terms.}
  \item{scale}{Scale the columns of X? Default is scale=TRUE.}
  \item{k}{Number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=10. }
  \item{plot.it}{Plot the cross-validation error as a function of \code{lambda}? Default is FALSE.}
   \item{verbose}{Print information on conflicting signs etc. Default is \code{verbose=FALSE}}
}
\details{

}
\value{
  \item{pcor}{estimated matrix of partial correlations.}
  \item{lambda.opt}{optimal value of \code{lambda} for each of the \code{ncol} regression models.}
}
\references{ N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks using Gaussian Graphical Models", preprint

  \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}

\seealso{\code{\link{ridge.cv}}}
\examples{
n<-50
p<-10
X<-matrix(rnorm(n*p),ncol=p)
pc<-ridge.net(X,k=5)
}
\keyword{multivariate}