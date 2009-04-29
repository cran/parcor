\name{performance.pcor}
\alias{performance.pcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality of estimated partial correlations}
\description{
  This function computed various performance measures of the estimated matrix of partial correlations.
}
\usage{
performance.pcor(inferred.pcor, true.pcor=NULL, fdr=TRUE, cutoff.ggm=0.8,verbose=FALSE,plot.it=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inferred.pcor}{matrix of estimated partial correlations}
  \item{true.pcor}{true matrix of partial correlations. Default is true.pcor=NULL}
  \item{fdr}{logical. If fdr=TRUE, the entries of \code{inferred.pcor} are tested for significance. Default is fdr=TRUE}
  \item{cutoff.ggm}{default cutoff for significant partial correlations. Default is cutoff.ggm=0.8}
  \item{verbose}{Print information on test results etc.. Default is \code{verbose=FALSE}}
  \item{plot.it}{Plot test results and ROC-curves. Default is \code{plot.it=FALSE}}
}
\details{This function computes a range of performance measures: The function always returns the number of selected edges, the binary matrix that encodes
the edges, the connectivity and the percentage of positive correlations. If \code{true.pcor} is provided, the function also returns
the power (= true positive rate), the false positive rate and the positive predictive value. For non-sparse estimates that involve
testing (i.e. \code{fdr=TRUE}) the function also returns the area under the curve, and a pair of vectors of false and true positive rates.
The latter can e.g. be used to plot a ROC-curve.
}

\value{
  \item{num.selected}{number of selected edges}
  \item{adj}{binary matrix that encodes the existence of an edge between two nodes.}
  \item{connectivity}{vector of length \code{ncol(inferred.pcor)}. Its ith entry indicated the number of nodes that are connected to the ith node.}
  \item{positive.cor}{percentage of positive partial correlations out of all selected edges.}
  \item{power}{power (if true.pcor is provided)}
  \item{ppv}{positive predictive value (if true.pcor is provided)}
  \item{tpr}{true positive rate (=power) (if true.pcor is provided) }
  \item{fpr}{true positive rate (=power) (if true.pcor is provided)}
  \item{auc}{area under the curve (if true.pcor is provided and \code{fdr=TRUE})}
  \item{TPR}{vector of true positive rates corresponding to varying cut-offs (if true.pcor is provided and \code{fdr=TRUE})}
  \item{FPR}{vector of false positive rates corresponding to varying cut-offs (if true.pcor is provided and \code{fdr=TRUE})}
}
\author{Juliane Schaefer, Nicole Kraemer}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
