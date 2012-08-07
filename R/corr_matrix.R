## Creating the correlation matrix based on hyperparameters 
## beta and a data matrix x
##
## May 8th, 2012

corr_matrix <- function(X,beta){
## Checking to make sure the data is a matrix, and sets it as one 
## if it is not
if (is.matrix(X) == FALSE){
	X = as.matrix(X)
}
d = ncol(X)
n = nrow(X)
## Checking the dimensions between the two inputs 
if (d != length(beta)){
	stop("The dimensions of beta and X do not match. \n")
}
dim(beta) = c(d,1)
R = matrix(0,n,n)
Points = expand.grid((1:n),(1:n))
Points = cbind(Points[,2],Points[,1])
Points = Points[Points[,2]>Points[,1],]
junk = X[Points[,2],] - X[Points[,1],]
junk = junk*junk
Beta = matrix(beta,nrow = (length(junk)/d),ncol=d,byrow=TRUE)
Rtemp=(10^Beta)*junk
Rtemp=rowSums(Rtemp)
R[Points]=Rtemp
R = R+t(R)
R=exp(-R)
return(R)
}

