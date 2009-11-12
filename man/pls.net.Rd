\name{pls.net}
\alias{pls.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Partial Correlations with Partial Least Squares}
\description{
 This function computes the matrix of partial correlations via an
  estimation of the corresponding regression models via Partial Least Squares.
}
\usage{
pls.net(X, scale = TRUE, k = 10, ncomp = 15,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables.}
  \item{scale}{Scale the columns of X? Default is scale=TRUE.}
  \item{k}{Number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=10. }
  \item{ncomp}{Maximal number of components. Default is 15.}
   \item{verbose}{Print information on conflicting signs etc. Default is \code{verbose=FALSE}}
}
\details{
 For each of the columns of \code{X}, a regression model based on
  Partial Least Squares is computed. The optimal model is determined via
  cross-validation.  The results of the regression models are
  transformed via the function \code{Beta2parcor}.
}
\value{
  \item{pcor}{estimated matrix of partial correlation coefficients.}
  \item{m}{optimal number of components for each of the \code{ncol(X)} regression models.}
}
\references{
 N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks using Gaussian Graphical Models", BMC Bioinformatics, in press}

\author{Nicole Kraemer}


\examples{
n<-20
p<-40
X<-matrix(rnorm(n*p),ncol=p)
pc<-pls.net(X,ncomp=10,k=5)
}
\keyword{multivariate}
