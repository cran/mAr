"mAr.eig" <-
function(A,C, ...)
{

# A has m rows and mp columns
m=dim(A)[1]
p=(dim(A)[2])/m

# Augmented coefficient matrix
At=matrix(nrow=m*p,ncol=m*p)
if (p==1) At=A else{
At[seq(1,m),seq(1,m*p)]=A
At[seq(m+1,m*p),seq(1,m*p-m)]=diag(1,(p-1)*m)
At[seq(m+1,m*p),seq(m*p-m+1,m*p)]=0
}

# Computation of eigenvalues and eigenvectors
l=eigen(At)$values
V=eigen(At)$vectors

# test stability of AR model
if ((any(Mod(l)>1))) warning("unstable AR model")


# fix phase of eigenvectors to satisfy normalisation conditions
a=matrix(nrow=1,ncol=dim(V)[2])
b=matrix(nrow=1,ncol=dim(V)[2])
St=matrix(nrow=dim(V)[2],ncol=dim(V)[2])
for (j in seq(1,dim(V)[2])){
a=Re(V[,j])
b=Im(V[,j])
ph=0.5*atan(2*sum(a*b)/( b %*% b - a %*% a))
na=sqrt(sum((cos(ph)*a - sin(ph)*b)^2))
nb=sqrt(sum((sin(ph)*a + sin(ph)*b)^2))
if(nb>na && ph < 0){ph=ph-pi/2}
if(nb>na && ph >0){ph=ph+pi/2}
St[,j]=V[,j] %*% exp(1i*ph)
}

# m-dimensional eigenvectors
S=St[seq(1+(p-1)*m,p*m),]

# transformed noise-covariance matrix
StInv=solve(St)[,seq(1,m)]
Ct=StInv %*% C %*% t(Conj(StInv))

# damping time and period for eigenmode i
tau=matrix(nrow=1,ncol=m*p)
per=matrix(nrow=1,ncol=m*p)
exctn=matrix(nrow=1,ncol=m*p)
for (i in seq(1,m*p)){
tau[i]=-2/log((abs(l[i]))^2)
a=Re(l[i])
b=Im(l[i])
if (identical(b,0)  && a>=0) {per[i]=Inf}
if (identical(b,0)  && a<0) {per[i]=2}
else {per[i]=2*pi/abs(atan2(b,a))}
exctn[i]=Re(Ct[i,i]/(1-(abs(l[i]))^2))
}

# relative dynamical importance of modes
exctn=exctn/sum(exctn)

return (list(dampTime=tau,period=per,excitations=exctn,eigv=S))

}
