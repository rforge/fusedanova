##' Fit a Fused ANOVA model
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total Variation) penality, 
##' ie. a sum of weighted \eqn{\ell_1}{l1}-norm on the difference of each coefficient.
##'
##' @param x matrix which rows represent individuals and columns 
##' independant variables. 
##'
##' @param class vector or factor giving the class of each individual.
##'
##' @param ... list of additional parameters to overwrite the defaults of the
##' fitting procedure. Include :
##' \itemize{%
##'
##' \item{\code{weights}: } {character; which type of weights is supposed to be used.
##' The supported weights are : \code{"default"}, \code{"laplace"}, \code{"gaussian"},
##' \code{"adaptive"}, \code{"naivettest"}, \code{"ttest"}, \code{"welch"} and \code{"personal"}. 
##' By default, its value is \code{"default"}.}
##'
##' \item{\code{gamma}: } {numeric; the \eqn{\gamma}{gamma} parameter needed for 
##' \code{"laplace"}, \code{"gaussian"} and \code{"adaptive"} weights. By default, 0.}
##' 
##' \item{\code{W}: } {numeric matrix; the matrix of weights needed if the \code{"personal"}
##' weights were selected. By default, a matrix with zero row and zero column.}
##' 	
##' \item{\code{standardize}: }{ logical; should each variable be standardized before the calculus ?
##' By default, \code{TRUE}.
##' }	
##'
##'	\item{\code{splits}: }{ integer; coding for forcing split or nosplit algorithms :
##' \itemize{%
##'
##' \item{\code{0} : }{Default value, let the programm decide 
##' which algorithm to use depending on the choosen \code{weights}.}
##' 
##' \item{\code{1} : }{Forces the algorithm not to take the splits into account.}
##'
##' \item{\code{2} : }{Forces the algorithm to take the splits into account even if 
##' the solution paths contains no split.}
##' }
##' Note : For the moment, only the no split algorithm has been coded. Please ensure that 
##' your weights choice leads to the no split algorithm or set \code{split} to \code{1}.
##' }
##'
##'	\item{\code{epsilon}: }{numeric; tolerance parameter for numeric calculations.
##' By default, \eqn{10^-10}{eps}. 
##' Note : this is currently not used.}
##' 
##' \item{\code{checkargs}: }{logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is \code{TRUE}. 
##' Automatically set to \code{FALSE} when a call is made
##' from cross-validation}		
##'	
##' \item{\code{lambdalist}: }{numeric vector; a set of \eqn{\lambda}{lambda} value for which
##' a prediction is asked. By default, a null vector. If the length of \code{lambdalist} 
##' is not \code{0}, the \code{fusedanova} class returnded by the \code{fusedanova} function will
##' have a not null attribute \code{prediction}. 	
##' }
##'
##' \item{\code{mc.cores}: } {integer; the number of cores to use. The default uses all
##' the cores available. }
##'
##' \item{\code{verbose}: } {boolean; should the code print out its progress.
##' By default, FALSE. }
##'
##' \item{\code{mxSplitSize}: } {integer; the maximum size for a group for being checked
## for eventual splitting. By default, 100. 
##' the cores available. }
##'
##' }
##'
##' @return an object with class \code{fusedanova}, see the
##' documentation page \code{\linkS4class{fusedanova}} for details.
##'
##' @note The optimized criterion is the following: 
##'
##' @seealso See also \code{\linkS4class{fusedanova}},
##' \code{\link{plot.fusedanova}} and \code{\link{crossval}}.
##' @name fusedanova
##' @rdname fusedanova
##' @keywords models, regression
##'
##' @examples
##' data(aves)
##' fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
##' plot(fa.laplace, labels=aves$order)
##' 
##' fa.ttest <- fusedanova(x=aves$weight, class=aves$family, weights="naivettest")
##' plot(fa.ttest, labels=aves$order)
##' 
##' fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
##' plot(fa.ada, labels=aves$order)
##' 
##' @export
fusedanova <- function(x,class,...) {

	## ===================================================
	## GETTING DEFAULT PARAMS
	user <- list(...)
	defs <- default.args.fa()
	args <- modifyList(defs, user)

	if (Sys.info()[['sysname']] == "Windows") {	
		args$mc.cores <- 1 # Windows does not support fork
	}
	
	## ===================================================
	## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
        if(!inherits(x, c("matrix", "Matrix")))
          x <- as.matrix(x)
	p <- ncol(x) # problem size
	n <- nrow(x) # sample size

	if (args$checkargs) {
                if(!inherits(x, c("matrix", "Matrix")))
			stop("x has to be of class 'matrix'.")
		if(any(is.na(x)))
			stop("NA value in x not allowed.")
		if(n != length(class))
			stop("x and y have not correct dimensions")
		if(length(unique(class))==1)
			stop("y has only one level.") 
		if(!(args$weights %in% possibleWeights))
			stop("Unknown weight parameter formulation. Aborting.") 
		# if (args$weights %in% varNeededWeights && sum(tabulate(class)==1))
			# stop("Weights requires variance. There is a class of length 1. Aborting.") 
		if (Sys.info()[['sysname']] == "Windows") {
			if(args$verbose){
				warning("\nWindows does not support fork, enforcing mc.cores to '1'.")	
			}
		}
	}
	
	if (length(args$lambdalist)!=0){ # sort by descending for effective prediction
		args$lambdalist = sort(args$lambdalist, decreasing = TRUE)
	}


	## ====================================================
	## CONVERSIONS OF x in a list and y in a numeric vector <================================= ?????????????????? should not be needed. toby : LAPPLY(1:ncol(x),function(k){ => ok for lapply 
	x = split(x, rep(1:p, each = n))
	if(!is.factor(class)){class=as.factor(class)}
	
	
	## ======================================================
	## NORMALIZATION TREATMENT
	if (args$standardize==TRUE){x=lapply(x,function(z){normalize(z,class)})} # centrer réduire

	## ====================================================
	## SPLIT OCCURENCE TESTING AND LAUNCH
	if (args$splits==0){ # if we let the programm choose, overwrite the value
		if (args$weights=="default"||(args$weights %in% c("laplace","gaussian","adaptive") && args$gamma>=0) ||length(unique(class))<3){
			args$splits = 1
		}else{
			args$split=2	
		}
	}
	
	if (args$split==1){
		algoType = "No Split"
		if(args$verbose){print("Path calculated with path without split")}
	}else{
		algoType = "With possible Splits"
		if(args$verbose){print("Path calculated with path with possible splits")}
	}
	
	result = mclapply(x,function(z) {calculatepath(z,class,args)},
          mc.cores= args$mc.cores, mc.preschedule=ifelse(p > 100,TRUE,FALSE)) # path algorithm 
        
	if (length(args$lambdalist) == 0){ # default parameter
		l = numeric(0)
		prediction =  list()
	}else{
		l =args$lambdalist
		prediction = lapply(result, function(z){list(table=z$prediction, order = z$order)})
	}

	result = lapply(result, function(z){list(table=z$table, order = z$order)})
	
	# small warning on last beta
	if (args$standardize==TRUE && abs(result[[1]]$table[1,1])>10^(-8)){
	warning("There may be some approximation errors (Beta(lambda_max) far from 0). You may want to lower the gamma if you are using one.")} 
 
	return(new("fusedanova",result=result, prediction=prediction, classes = class, weights=args$weights, lambdalist=l,algorithm=algoType))

}


#########################################
# Calculate path with or without splits #
#########################################
# calculate the path for one gene
# x the vector of data, group the vector group belonging, args all the rest
calculatepath <- function(x,group,args){

	ngroup = tabulate(group)# vector of number by group
	xm =  rowsum(x,group)/ngroup
	xv = rep(0,length(xm))
	if (args$weights %in% varNeededWeights){
		# var needed if weights are of welch or ttest type
		xv = ngroup/(ngroup-1)*(rowsum(x^2,group)/ngroup - xm^2)
	}
 
	ordre = order(xm)
	xm = xm[ordre] # sort from the smallest beta to the highest
	ngroup = ngroup[ordre]
	xv =xv[ordre]
	
	if (args$splits==1){
		res  <- .Call("noSplit",R_x=xm,R_xv=xv,R_ngroup=ngroup, R_args=args, PACKAGE="fusedanova")
	}else{
		res  <- .Call("withSplit",R_x=xm,R_xv=xv,R_ngroup=ngroup, R_args=args, PACKAGE="fusedanova")
	}

	return(list(table=res$res, prediction = res$pred, order=ordre))
  
} 

#############################
# normalization of a vector
#############################
normalize <- function(x,group){
	Pooled = (nlevels(group)!= length(x))   
	if (Pooled ==FALSE){  
		res = (x-mean(x))/as.numeric((sqrt(var(x))))
	}else{
		n = length(x)
		m  = mean(x)
        ngroup = tabulate(group)
		s =  1/(ngroup-1)*(rowsum(x^2,group) - (1/ngroup)*(rowsum(x,group))^2)
		s[which(ngroup==1)]=0
		if (sum(s)==0){
			res = (x-mean(x))/as.numeric((sqrt(var(x))))
		}else{
			s = sqrt(sum(s*(ngroup-1))/(n-length(ngroup)))
			res = (x-m)/s
		}
	}
	return(res)
}

#############################
# default args
#############################
default.args.fa <- function() {
	return(list(
		weights = "default",
		W = matrix(nrow=0, ncol=0),
		gamma = 0 ,
		standardize =TRUE,
		splits = 0,
		epsilon =10^-10,
		checkargs = TRUE,
		lambdalist = numeric(0),
		mc.cores = detectCores(),
		verbose = FALSE,
		mxSplitSize = 100
    ))
}

# constant list of treated weights
possibleWeights = c("default","laplace","gaussian","adaptive","naivettest","ttest","welch","personal")
varNeededWeights = c("welch", "ttest")
