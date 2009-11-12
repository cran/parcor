\name{Beta2parcor}
\alias{Beta2parcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of partial correlation coefficients}
\description{
  This function computes the matrix of partial correlation coefficients
  based on the results of the corresponding regression models.
}
\usage{
Beta2parcor(Beta,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Beta}{matrix of regression coefficients}
  \item{verbose}{print information on conflicting signs etc. Default is \code{verbose=FALSE}.}
}
\details{A well-known result (Whittaker, 1990) shows that the matrix of
  partial correlation coefficients can be estimated by computing a least
  squares regression model for each variable. If there are more
  variables than observations, the least squares problem is ill-posed
  and needs regularization. The matrix \code{Beta} stores the regression
  coefficients of any user-defined regression method. The function
  \code{Beta2parcor} computes the
  corresponding matrix of partial correlations.
}
\value{
  matrix of partial correlation coefficients
}
\references{J. Whittaker (1990)  "Graphical models in applied
  multivariate statistics", Wiley, New York.
  
 N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", BMC Bioinformatics, in press
}
\author{Nicole Kraemer}
\note{This is an internal function.

}
\seealso{\code{\link{ridge.net}}, \code{\link{adalasso.net}},\code{\link{pls.net}}}
\examples{
# this is an internal function and should not be called by the user
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
