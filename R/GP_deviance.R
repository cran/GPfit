## Calculates the deviance values 
##
## May 8th, 2012

GP_deviance <- function(beta,X,Y,nug_thres=20){
if (is.matrix(X) == FALSE){
	X = as.matrix(X)
}
n = nrow(X);
d = ncol(X);
if (n != length(Y)){
	stop("The dimensions of X and Y do not match. \n")
}
## Checking the dimensions of the parameters
if (d != length(beta)){
	stop("The dimensions of beta and X do not match \n")
}
if (n != length(Y)){
	stop("The dimensions of X and Y do not match \n")
}

if (nug_thres < 10){
	warning("nug_thres outside of the normal range of [10, 25]")
}
if (nug_thres > 25){
	warning("nug_thres outside of the normal range of [10, 25]")
}

One = rep(1,n);
dim(Y) = c(n,1);

dim(beta) = c(1,d);

R = corr_matrix(X,beta);
temp = eigen(R,symmetric = TRUE, only.values = TRUE);
eig_val = temp$values;

condnum = kappa(R,triangular = TRUE);
max_eigval = eig_val[1]
delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*
		(exp(nug_thres)-1))))

LO = diag(n);
Sig = R + delta*LO;
L = chol(Sig);

Sig_invOne = solve(L,solve(t(L),One));
Sig_invY = solve(L,solve(t(L),Y));
Sig_invLp = solve(L,solve(t(L),LO));

mu_hat = solve(t(One)%*%Sig_invOne, t(One)%*%Sig_invY);
Sig_invb = solve(L,solve(t(L),(Y-One%*%mu_hat)));

MM = max(max(Sig_invLp));
Sig_invLpp = Sig_invLp/MM;

part1_1 = 2*log(abs(det(L))) 
part1 = -n * log(MM) - log(abs(det(Sig_invLpp)))

if(is.finite(part1_1)){ part1 = part1_1;}

part2 = n * log(t(Y - One %*% mu_hat) %*% Sig_invb)
devval = part1 + part2

devval = devval[1,1];
if (is.finite(devval) == FALSE)
{
	stop('Infinite values of the Deviance Function, 
		unable to find optimum parameters \n')
}
return(devval)
}
