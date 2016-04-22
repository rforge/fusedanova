context("Consistency of the fused anova solution path with the flsa package")

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

test_that("Consistency between 'fusedanova' and 'flsa' packages", {

	require(flsa)

	# for the moment using random data
	splits = 1
	xm = as.matrix(runif(10,min=0,max=20))
	ng = rep(1,10)
	sigma = 1
	Y <- Sim(xm,ng,sigma)
	class= Y[,1]
	Y = as.matrix(Y[,2])
	
	connList = vector("list",10)
	for (i in 1:10){
		connList[[i]] = (0:9)[-i]
	} # in fusedanova the graph is always cliquish
	
	fa = fusedanova(x=Y,class=class,standardize=FALSE,weights = "default")
	falambda = sort(unique(fa@result[[1]]$table$lambda)) 
	
	fl <- flsa(Y, connListObj=connList)
	fllambda = sort(unique(fl$BeginLambda))
	
    expect_that(falambda, is_equivalent_to(fllambda))
	# ie. same lambda list 
	# fused anova should be equivalent to flsa (with clique graph) for 
	# default weights and innitial groups of size 1.
	
	
	# test of consistency with personal W in 0 or 1 to mimmic a graph and make it split.
	# if everything go well, 5 should split with 6,7,8,9
	connList = vector("list",9)
	connList[[1]] = as.integer((1:8))
	connList[[2]] = as.integer(c(0,(2:3),(5:8)))
	connList[[3]] = as.integer(c((0:1),3,(5:8)))
	connList[[4]] = as.integer(c((0:2),(5:8)))
	connList[[5]] = as.integer(c(0,8)) #<--------------- here is the tricky part
	connList[[6]] = as.integer((0:8)[-c(5,6)])
	connList[[7]] = as.integer((0:8)[-c(5,7)])
	connList[[8]] = as.integer((0:8)[-c(5,8)])
	connList[[9]] = as.integer((0:8)[-9])
	
	Y = c(1000,1001,1002,1003,10000,11000,11001,11002,11003)
	
	W = matrix(1,9,9)
	diag(W) =0
	W[5,2:8] = 0
	W[2:8,5]=0
	
	fl <- flsa(Y, connListObj=connList)
	fllambda = sort(unique(as.vector(fl$BeginLambda)))
	
	
	Y=as.matrix(Y)
	class= 1:9
	fa = fusedanova(x=Y,class=class,standardize=FALSE,weights = "personal",W = W)
	falambda = sort(unique(fa@result[[1]]$table$lambda)) 
	
	expect_that(falambda, is_equivalent_to(fllambda))
	
})
