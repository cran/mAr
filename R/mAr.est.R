"mAr.est" <-
function(x,pmin,pmax, ...)
{

x=as.matrix(x)
n=dim(x)[1]
m=dim(x)[2]

if(pmin >= pmax)
	stop("pmin must be < pmax")

if(round(pmin) != pmin || round(pmax) != pmax)
	stop("order model must be an integer")

ne=n-pmax	# length of time vector
npmax=m*pmax+1

if(ne <= npmax)
	stop("time series too short")


#
# Model order selection
#


# Model of order pmax

K=matrix(nrow=ne,ncol=npmax+m)
K[,1]=rep(1,ne)
for(j in 1:pmax) { K[,seq(2+m*(j-1),1+m*j)]=data.matrix(x[seq(pmax-j+1,n-j),1:m]) }
K[,seq(npmax+1,npmax+m)]=data.matrix(x[seq(pmax+1,n),1:m])

q=ncol(K)
delta=(q^2+q+1)*(.Machine$double.eps) # regularisation parameter
scale=sqrt(delta)*sqrt(apply(K^2,2,sum))

R=qr.R(qr((rbind(K,diag(scale)))),complete=TRUE) # QR factorisation of data matrix K

R22=R[seq(npmax+1,npmax+m),seq(npmax+1,npmax+m)]

logdp=c(pmin:pmax)
logdp[pmax]=2*log(abs(prod(diag(R22)))) # log determinant of residual cross-product matrix

sbc=c(pmin:pmax)
sbc[pmax]=logdp[pmax]/m - log(ne)*(ne-npmax)/ne


# Downdating - approximate order selection criteria for models of order pmax-1:pmin

p=seq((pmax),pmin,-1)
np=c(pmin:pmax)
np[p]=m*p+1

invdp= solve(R22) %*% t(solve(R22))

for(i in seq(pmax-1,pmin,-1)){
Rp=R[seq(np[i]+1,np[i]+m),seq(npmax+1,npmax+m)]
L=chol(t(diag(1,m,m) + Rp %*% invdp %*% t(Rp)))
N=solve(L) %*% Rp %*% invdp
invdp=invdp - (t(N) %*% N)
logdp[i]= logdp[i+1] + 2*log(abs(prod(diag(L))))
sbc[i]=logdp[i]/m - log(ne)*(ne-np[i])/ne # Schwartz's Bayesian Criterion
}

# selected optimal order
popt=pmin+(which.min(sbc)-1)
npopt=m*popt+1


#
# Parameters estimation for optimal order model
#

Ropt11=R[seq(1,npopt),seq(1,npopt)]
Ropt12=R[seq(1,npopt),seq(npmax+1,npmax+m)]
Ropt22=R[seq(npopt+1,npmax+m),seq(npmax+1,npmax+m)]

Ropt11[,1]=Ropt11[,1]*max(scale[2:(npmax+m)])/scale[1] # re-scaling R11's first column to improve condition

B=t(solve(Ropt11) %*% Ropt12) # Estimated augmented parameter matrix

w=B[,1]*max(scale[2:(npmax+m)])/scale[1] # intercept vector
A=B[,2:npopt] # coefficient matrix

C=(t(Ropt22) %*% Ropt22)/(ne-npopt) # bias-corrected covariance matrix estimate

# Residuals

t=seq(1,(n-popt))
res=matrix(nrow=(n-popt),ncol=m)
res[t,seq(1,m)]=x[t+popt,]-(rep(1,n-popt) %*% t(w))
for (j in seq(1,popt)) {res[t,seq(1,m)]=res[t,seq(1,m)]-(x[(t-j+popt),] %*% t(A[, seq(m*j-m+1,j*m)]))}


return (list(pHat=popt,SBC=sbc[pmin:pmax],wHat=w,AHat=A,CHat=C,res=res))


}
