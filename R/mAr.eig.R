"mAr.eig" <-
function(A,C, ...)
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

# Computation of eigenvalues and eigenvectors
V=eigen(At)
l=V$values

# test stability of AR model
if ((any(Mod(l)>1))) warning("unstable AR model")

# damping time and period for eigenmode i
tau=matrix(nrow=1,ncol=m*p)
per=matrix(nrow=1,ncol=m*p)
for (i in seq(1,m*p)){
tau[i]=-2/log((abs(l[i]))^2)
a=Re(l[i])
b=Im(l[i])
if (identical(b,0)  && a>=0) {per[i]=Inf}
if (identical(b,0)  && a<0) {per[i]=2}
else {per[i]=2*pi/abs(atan2(b,a))}
}

return (list(dampTime=tau,period=per))

}
