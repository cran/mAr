\name{mAr.eig}
\alias{mAr.eig}
\title{ Eigendecomposition of m-variate AR(p) model}
\description{
  Computation of the oscillation periods and damping rates of a system modeled by a m-variate AR(p) model from eigen-decomposition of the
estimated matrix of autoregressive coefficients
}
\usage{
mAr.eig(A, C, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ matrix of estimated autoregression coefficients}
  \item{C}{ noise covariance matrix}
  \item{\dots}{additional arguments for specific methods}
}

\value{
  A list with components: 
  \item{dampTime }{damping times associated with each eigenmode}
  \item{period }{periods associated with each eigenmode}
}

\references{
Neumaier, A. and Schneider, T. (2001), Estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Transactions on Mathematical Software, 27, 1, 27-57. \cr\cr
Schneider, T. and Neumaier, A. (2001), A Matlab package fo the estimation of parameters and eigenmodes of multivariate autoregressive models, 27, 1, 58-65. \cr\cr
}

\author{ S. M. Barbosa }

\examples{
data(pinkham)
y=mAr.est(pinkham,2,5)
mAr.eig(y$AHat,y$CHat)
}
\keyword{ multivariate}

