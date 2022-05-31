\name{mAr.eig}
\alias{mAr.eig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigendecomposition of m-variate AR(p) model }
\description{
 Eigen-decomposition of the estimated matrix of autoregressive coefficients from an m-variate AR(p) model
}
\usage{
mAr.eig(A, C = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ matrix of estimated autoregression coefficients }
  \item{C}{ noise covariance matrix }
  \item{\dots}{ additional arguments for specific methods}
}

\value{
  A list with components: 
  \item{modes }{periods and damping times associated to each eigenmode}
  \item{eigv}{m*p m-dimensional eigenvectors}
}
\references{ 
Barbosa S.M., Silva M.E., Fernandes M.J. (2006), Multivariate autoregressive modelling of sea level time series from TOPEX/Poseidon satellite altimetry. Nonlinear Processes in Geophysics, 13, 177-184.
\cr\cr
Neumaier, A. and Schneider, T. (2001), Estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Transactions on Mathematical Software, 27, 1, 27-57. \cr\cr
Schneider, T. and Neumaier, A. (2001), A Matlab package fo the estimation of parameters and eigenmodes of multivariate autoregressive models, 27, 1, 58-65. \cr\cr
}

\author{S. M. Barbosa }



\examples{
data(pinkham)
y=mAr.est(pinkham,2,5)
mAr.eig(y$AHat,y$CHat)

}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
