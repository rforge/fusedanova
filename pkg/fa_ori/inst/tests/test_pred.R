context("Consistency of the fused anova predictions")

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

test_that("Consistency between 'fusedanova' predictions", {

	# using random data
	# 1 and 2 -> no split algo, 3 and 4 -> splits are possible.
	
	xm = as.matrix(runif(10,min=0,max=20))
	ng = sample(1:100,size=10,replace=TRUE)
	sigma = 1
	Y <- Sim(xm,ng,sigma)
	class= Y[,1]
	Y = as.matrix(Y[,2])
	test1 = fusedanova(x=Y,class=class,standardize =TRUE,weights = "default") # classic
	test2 = fusedanova(x=Y,class=class,standardize =TRUE,weights = "default",lambdalist=seq(0,0.5,length.out=50)) # with lambdalist
	# prediction
	p1 = predict(test1,lambda=seq(0,0.5,length.out=50))
	p2 = predict(test2)
	
	test3 = fusedanova(x=Y,class=class,standardize =TRUE,weights = "naivettest") # classic
	test4 = fusedanova(x=Y,class=class,standardize =TRUE,weights = "naivettest",lambdalist=seq(0,5,length.out=50)) # with lambdalist
	# prediction
	p3 = predict(test3,lambda=seq(0,5,length.out=50))
	p4 = predict(test4)
	
	expect_that(p1["beta"], is_equivalent_to(p2["beta"]))
	expect_that(p3["beta"], is_equivalent_to(p4["beta"]))
	plot(test3)
	# p3[abs(p3["beta"]-p4["beta"])>10^-6,]
	
})

