"mAr.sim" <-
function(w,A,C,N, ...)
{

# A has m rows and mp columns
m=dim(A)[1]
p=(dim(A)[2])/m

# Augmented coefficient matrix
At=matrix(nrow=m*p,ncol=m*p)
if (p==1) {At=A}
else{
At[seq(1,m),seq(1,m*p)]=A
At[seq(m+1,m*p),seq(1,m*p-m)]=diag(1,(p-1)*m)
At[seq(m+1,m*p),seq(m*p-m+1,m*p)]=0
}

# test stability of AR model
l=(eigen(At,only.values = TRUE))$values
if ((any(Mod(l)>1))) warning("unstable AR model")

# discard first nd simulated values
nd=1000

U=chol(C) #upper triangular matrix

# generate independent gaussian noise vectors with mean 0 and covariance matrix C
# uses mvrnorm in MASS to simulate from a multivariate normal distribution
require(MASS)
noisevec=mvrnorm(nd+N, rep(0,m),C)
matw=rep(1,nd+N)%*%t(w)
vec=noisevec+matw

# p initial values equal to the process mean
if (any(w!=0)){
B=diag(1,m)
for(j in seq(1,p)) {B=B-A[,seq(m*j-m+1,j*m)]}
mproc=as.vector(solve(B) %*% w) # process mean
xi=(matrix(1,nrow=p,ncol=1)) %*% mproc
} else { xi=matrix(nrow=p,ncol=m)
xi[,]=0
}
# initialisation of state vector
u=matrix(nrow=p+nd+N,ncol=m)
u[seq(1,p),seq(1,m)]=xi
u[seq(p+1,p+nd+N),seq(1,m)]=0

# simulate nd+N observations of multivariate AR(p) process
Atr=t(A)
x=matrix(ncol=m,nrow=p)

for(k in seq(p+1,nd+N+p)){
for(j in seq(1,p)){
x[j,]=u[k-j,]%*%Atr[seq(m*j-m+1,m*j),]}
u[k,]=as.matrix(apply(x,2,sum)+ vec[k-p,])
}

# return N simulated observations
v=u[seq(nd+p+1,nd+p+N),]
simulated=data.frame(v[,seq(1,m)])

return (simulated)

}
