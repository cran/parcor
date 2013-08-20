\name{sym2vec}
\alias{sym2vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform symmetric matrix to vector}
\description{
  This function transforms the upper triangle (without diagonal) of a matrix A into a vector.
  
}
\usage{
sym2vec(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{quadratic matrix of size \code{p x p}}
}
\details{The entries of the matrix \code{A} are stacked line-by-line.
}
\value{
  vector of length p(p-1)/2
}

\author{Nicole Kraemer}
\note{This is an internal function.

}

\examples{
p<-7
A<-matrix(rnorm(p*p),ncol=p)
v<-sym2vec(A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
