\name{parcor-package}
\alias{parcor-package}
\alias{parcor}
\docType{package}
\title{Parcor: Estimation of partial correlations based on regularized regression.
}
\description{
This package contains different methods to estimate the matrix of
partial correlations based on a (n x p) matrix X. For p>n, the matrix of
partial correlations can be estimated based on p least-squares
regression fits. However, for p<n, theses least-squares problems are
ill-posed and need to be regularized. This package contains four different
regularized regression techniques for the estimation of the partial
correlations: lasso, adaptive lasso, ridge regression, and partial least squares. 

}
\details{
\tabular{ll}{
Package: \tab parcor\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-04-22\cr
License: \tab GPL2 or newer\cr
LazyLoad: \tab yes\cr
}

}
\author{Nicole Kraemer, Juliane Schaefer

Maintainer: Nicole Kraemer <nkraemer@cs.tu-berlin.de>
}
\references{
N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", preprint

  \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}
}
\keyword{ package }
