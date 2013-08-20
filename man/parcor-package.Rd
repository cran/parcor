\name{parcor-package}
\alias{parcor-package}
\alias{parcor}
\docType{package}
\title{Parcor: Estimation of partial correlations based on regularized regression.
}
\description{
This package contains different methods to estimate the matrix of
partial correlations based on a (n x p) matrix X of observation. For low-dimensional settings (p>n), the matrix of
partial correlations can be estimated based on p least-squares
regression fits. However, in high-dimensional scenarios (p<n), theses least-squares problems are
ill-posed and need to be regularized. This package contains four different
regularized regression techniques for the estimation of the partial
correlations: lasso, adaptive lasso, ridge regression, and Partial Least Squares. In addition, the package provides
model selection for lasso, adaptive lasso and Ridge regression based on cross-validation.

}
\details{
\tabular{ll}{
Package: \tab parcor\cr
Type: \tab Package\cr
Version: \tab 0.2-4\cr
Date: \tab 2013-08-20\cr
License: \tab GPL2 or newer\cr
LazyLoad: \tab yes\cr
}

}
\author{Nicole Kraemer, Juliane Schaefer

Maintainer: Nicole Kraemer <nicole.kraemer@tum.de>
}
\references{
N. Kraemer, J. Schaefer, A.-L. Boulesteix (2009) "Regularized Estimation of
  Large-Scale Gene Regulatory Networks with Gaussian Graphical Models", BMC Bioinformatics, 10:384

\url{http://www.biomedcentral.com/1471-2105/10/384/}
}
\keyword{ package }
