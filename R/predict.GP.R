## Prediction for gaussian process models
##
## May 8th, 2012

predict.GP <- function(object,xnew = object$X,...){
if (is.GP(object) == FALSE){
	stop("The object in question is not of class \"GP\" \n")
}
if (is.matrix(xnew) == FALSE){
	xnew = as.matrix(xnew)
}

X = object$X
Y = object$Y
n = nrow(X)
d = ncol(X)
if (d != ncol(xnew)){
	stop("The training and prediction data sets are of 
	different dimensions. \n")
}

beta = object$beta;
sig2 = object$sig2;
delta = object$delta;

dim(beta) = c(d,1);
R = corr_matrix(X,beta);

One = rep(1,n);
LO = diag(n);
Sig = R + delta*LO;
L = chol(Sig);

Sig_invOne = solve(L,solve(t(L),One));
Sig_invY = solve(L,solve(t(L),Y));
Sig_invLp = solve(L,solve(t(L),t(LO)));

nnew = nrow(xnew);
Y_hat = rep(0,nnew);
MSE = rep(0,nnew);
for(kk in 1:nnew)
{
	xn = matrix(xnew[kk,],nrow=1);
	r=exp(-((X-as.matrix(rep(1,n))%*%(xn))^2)%*%(10^beta));
	yhat_x = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invY;
	Y_hat[kk] = yhat_x;
	
	Sig_invr = solve(L,solve(t(L),r));
	
	cp_delta_r = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invr;
	
	cp_delta_Lp = (((1-t(r)%*%Sig_invOne)/(t(One)%*%Sig_invOne))%*%t(One)+t(r))%*%Sig_invLp;
	mse = sig2*(1-2*cp_delta_r+cp_delta_Lp%*%R%*%t(cp_delta_Lp));
	MSE[kk] = mse*(mse>0);
}

prediction = NULL

names = c()
for (i in 1:d)
{
	names[i] = paste("xnew.",i, sep="")
}
names[d+1] = "Y_hat"
names[d+2] = "MSE"

full_pred = cbind(xnew, Y_hat, MSE)
colnames(full_pred) = names
prediction$Y_hat = Y_hat
prediction$MSE = MSE
prediction$complete_data = full_pred
return(prediction)
}
