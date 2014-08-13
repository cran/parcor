\name{lm.ridge.univariate}
\alias{lm.ridge.univariate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ridge Regression for a single predictor variable}
\description{
  This function computes the  ridge regression coefficients for a mdoel with a single predictor variable.
}
\usage{
lm.ridge.univariate(x,y,lambda=0,scale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of predictor observations.}
  \item{y}{vector of responses. The length of y must equal the length of x}
  \item{lambda}{vector of penalty terms. Default is lambda=0}
  \item{scale}{Scale x and y? Default is scale=TRUE.}
}

\details{This function replaces the R function \code{lm.ridge} if only one predictor variable is used, as the latter function produces an error message in this case.}

\value{A matrix with two columns and \code{length(lambda)} rows. In each row, the intercept (first colum) and the regression coefficient (second column) is stored.}

\author{Nicole Kraemer}

\seealso{\code{\link{ridge.net}},\code{\link{ridge.cv}} }
\examples{
n<-100 # number of observations
x<-rnorm(100)
y<-rnorm(n)
ridge.object<-lm.ridge.univariate(x,y,lambda=1:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univariate}
