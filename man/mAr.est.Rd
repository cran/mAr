\name{mAr.est}
\alias{mAr.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of multivariate AR(p) model }
\description{
Stepwise least-squares estimation of a multivariate AR(p) model based on the algorithm of Neumaier and Schneider (2001).
}
\usage{
mAr.est(x, p, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of multivariate time series}
  \item{p}{ model order}
  \item{\dots}{ additional arguments for specific methods }
}
\details{
   Fits by stepwise least squares an m-variate AR(p) model given by
  \deqn{X[t]=w + A1 X[t-1] +...+ Ap X[t-p] +e[t]}
  where \cr
  X[t]=[X1(t)...Xm(t)]' is a vector of length m \cr
  w is a m-length vector of intercept terms \cr
  A=[A1 ... Ap] is a mp x m matrix of autoregressive coefficients \cr
  e(t) is a m-length uncorrelated noise vector with mean 0 and m x m covariance matrix C
}
\value{
  A list with components:
  \item{SBC }{Schwartz Bayesian Criterion}
  \item{wHat }{vector of intercept terms}
  \item{AHat }{matrix of estimated autoregression coefficients for the fitted model}
  \item{CHat }{noise covariance matrix}
  \item{resid }{residuals from the fitted model}
}
\references{
Neumaier, A. and Schneider, T. (2001), Estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Transactions on Mathematical Software, 27, 1, 27-57. \cr\cr
Schneider, T. and Neumaier, A. (2001), A Matlab package fo the estimation of parameters and eigenmodes of multivariate autoregressive models, 27, 1, 58-65. \cr\cr
Lutkepohl, H. (1993), Introduction to Multiple Time Series Analysis. Springer-Verlag, Berlin.
}

\author{ S. M. Barbosa  }

%\seealso{
%        \code{\link[ts]{ar}}
%}

\examples{
data(pinkham)
y=mAr.est(pinkham,2,5)
}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
