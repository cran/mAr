\name{mAr.pca}
\alias{mAr.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate autoregressive analysis in PCA space}
\description{
Estimation of m-variate AR(p) model in reduced PCA space (for dimensionality reduction) and eigen-decomposition of augmented coefficient matrix 
}
\usage{
mAr.pca(x, p, k = dim(x)[2], ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of multivariate time series}
  \item{p}{ model order }
  \item{k}{ number of principal components to retain }
  \item{\dots}{ additional arguments for specific methods }
}

\value{
  A list with components:
  \item{p}{model order}      
  \item{SBC }{Schwartz Bayesian Criterion}
  \item{fraction.variance}{fraction of variance explained by the retained components}
  \item{resid}{residuals from the fitted model}
  \item{eigv}{m*p m-dimensional eigenvectors}
  \item{modes}{periods and damping times associated to each eigenmode}
}
\references{ Neumaier, A. and Schneider, T. (2001), Estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Transactions on Mathematical Software, 27, 1, 27-57. \cr\cr}

\author{ S. M. Barbosa}

\seealso{
        \code{\link[mAr]{mAr.est}}
}
\examples{

data(sparrows)
A=mAr.est(sparrows,1)$AHat
mAr.eig(A)$modes
mAr.pca(sparrows,1,k=4)$modes

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS

