context("Consistency of the fused anova solution path with the clusterpath")

test_that("Consistency between 'fusedanova' and 'clusterpath' packages", {

	require(clusterpath)

	# iris data :
	breakpoints  <- clusterpath.l1.id(iris[,1:4])
	clustlambda = unique(subset(breakpoints,col == "Sepal.Length" )$lambda)
	
	fa = fusedanova(x=as.matrix(iris[,1:4]),class=(1:length(iris[,1])),standardize =FALSE,weights = "default",checkargs = FALSE)
	falambda = unique(fa@result[[1]]$table$lambda) 
	
    expect_that(clustlambda, is_equivalent_to(falambda))
	# ie. same lambda list for the first variable.
	# fused anova should be equivalent to clusterpath for 
	# default weights and innitial groups of size 1.
	
})
