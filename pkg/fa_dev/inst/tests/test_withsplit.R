context("Consistency between 'fusedanova' with and without splits")

Sim <- function(xm,ngroup,sigma){
# xm mean of the group, ngroup number of person in the group, sigma vector or variances
# nb of lines = nb of groups
	p = ncol(xm)
	n = sum(ngroup)  
	Y = matrix(0,n,p+1)
	if(!inherits(sigma, c("matrix", "Matrix"))){
		sigma = matrix(sigma,nrow=nrow(xm), ncol=p) 
	}
	l=1
	for (i in 1:nrow(xm)){
		for(j in 1:ngroup[i]){ 
			Y[l,1] = i # numb of the group
			for (k in 1:p){	
					Y[l,k+1] = xm[i,k]+rnorm(1,0,sigma)
			}
			l=l+1
		}
	}
	return(Y)
}

test_that("Consistency between 'fusedanova' with and without splits", {

	# for the moment using random data
	weights = "default" 
	xm = as.matrix(runif(10,min=0,max=20))
	ng = sample(1:100,size=10,replace=TRUE)
	sigma = 1
	Y <- Sim(xm,ng,sigma)
	class= Y[,1]
	Y = as.matrix(Y[,2])
	
	fa = fusedanova(x=Y,class=class,weights = weights,checkargs = FALSE, splits=1)
	lambda1 = sort(unique(fa@result[[1]]$table$lambda)) 
	
	fa =  fusedanova(x=Y, class=class, weights = weights,checkargs = FALSE, splits=2, mxSplitSize=1)
	lambda2 = sort(unique(fa@result[[1]]$table$lambda)) 
	
	fa = fusedanova(x=Y, class=class, weights = weights,checkargs = FALSE, splits=2)
	lambda3 = sort(unique(fa@result[[1]]$table$lambda)) 

	expect_that(lambda1, is_equivalent_to(lambda2))
	# ie. same lambda list for the first variable. if we never check the mx flow
	
	expect_that(lambda3, is_equivalent_to(lambda2))
	# ie. same lambda list for the first variable. if we check the maxflow
	
	
})
