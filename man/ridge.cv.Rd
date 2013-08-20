\name{ridge.cv}
\alias{ridge.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ridge Regression.}
\description{
  This function computes the optimal ridge regression model based on cross-validation.
}
\usage{
ridge.cv(X, y, lambda, scale = TRUE, k = 10, plot.it = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input observations. The rows of \code{X} contain the
    samples, the columns of \code{X} contain the observed variables}
  \item{y}{vector of responses. The length of y must equal the number of
  rows of X}
  \item{lambda}{Vector of penalty terms.}
  \item{scale}{Scale the columns of X? Default is scale=TRUE.}
  \item{k}{Number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=10. }
  \item{plot.it}{Plot the cross-validation error as a function of \code{lambda}? Default is FALSE.}
}


\value{
  \item{intercept}{cross-validation optimal intercept}
  \item{coefficients}{cross-validation optimal regression coefficients}
  \item{lambda.opt}{optimal value of \code{lambda}.}
}
\author{Nicole Kraemer}

\seealso{\code{\link{ridge.net}}}
\examples{
n<-100 # number of observations
p<-60 # number of variables
X<-matrix(rnorm(n*p),ncol=p) 
y<-rnorm(n)
ridge.object<-ridge.cv(X,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
