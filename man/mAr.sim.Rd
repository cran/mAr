\name{mAr.sim}
\alias{mAr.sim}
\title{ Simulation from a multivariate AR(p) model }
\description{
  Simulation from an m-variate AR(p) model
}
\usage{
mAr.sim(w, A, C, N, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{ vector of intercept terms }
  \item{A}{ matrix of AR coefficients }
  \item{C}{ noise covariance matrix }
  \item{N}{ length of output time series }
  \item{\dots}{ additional arguments}
}
\details{
  Simulation from an m-variate AR(p) model given by 
  \deqn{X[t]=w + A1 X[t-1] +...+ Ap X[t-p] +e[t]}
  where \cr
  X[t]=[X1(t)...Xm(t)]' is a vector of length m \cr
  w is a m-length vector of intercept terms \cr
  A=[A1 ... Ap] is a m x mp matrix of autoregressive coefficients \cr
  e(t) is a m-length uncorrelated noise vector with mean 0 and m x m covariance matrix C
}
\value{
  returns a list containg the N simulated observations for each of the m time series
}

\references{Neumaier, A. and Schneider, T. (2001), Estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Transactions on Mathematical Software, 27, 1, 27-57. \cr\cr
Schneider, T. and Neumaier, A. (2001), A Matlab package fo the estimation of parameters and eigenmodes of multivariate autoregressive models, 27, 1, 58-65. \cr\cr
Lutkepohl, H. (1993), Introduction to Multiple Time Series Analysis. Springer-Verlag, Berlin.
 } 

\author{ S. M. Barbosa}


\examples{
w=c(0.25,0.1)
C=rbind(c(1,0.5),c(0.5,1.5))
A=rbind(c(0.4,1.2,0.35,-0.3),c(0.3,0.7,-0.4,-0.5))
x=mAr.sim(w,A,C,N=300)
}

\keyword{ multivariate }

