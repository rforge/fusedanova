##' Cross-validation function for fusedanova method.
##'
##' Function that computes K-fold cross-validated error of a
##' \code{fusedanova} fit. Possibility to perform V times the whole
##' K-fold CV to average the error and reduce the reduce the CV
##' standard error around the minimum.
##'
##' @param x matrix whose rows represent individuals and columns
##' independant variables.
##'
##' @param class vector or factor giving the class of each individual.
##'
##' @param K integer indicating the number of folds. Default is 10.
##'
##' @param V integer indicating the number of times the K folds CV
##' will be averaged. Default is 1 (no averaging).
##'
##' @param verbose boolean for verbose mode. Default is \code{TRUE}.
##'
##' @param folds list of \code{K} vectors that describes the folds to
##' use for the cross-validation. By default, the folds are randomly
##' sampled with the specified K. The same folds are used for each
##' values of \code{lambdalist}.
##'
##' @param lambdalist list of \code{lambda} penalty parameters used
##' in the cross validation. By default, lambdalist is NULL and calculated
##' using the maximum \code{lambda} and the parameter \code{nlambda}.
##'
##' @param ... list of additional parameters to overwrite the defaults of the
##' fitting procedure. See the corresponding documentation (\code{\link{fusedanova}}).
##' Also include :
##' \itemize{%
##'
##' \item{\code{nlambda}: } {integer; the length
##' of the \code{lambdalist} vector, by default 100}
##'
##' \item{\code{log.scale}: } {boolean; should a logarithmic scale be used
##' during the creation of \code{lambdalist}. By default, FALSE.
##' }
##'
##' \item{\code{min.ratio} : }{ numeric parameter setting the smallest value of lambdalist
##' in the \code{log.scale} case with the formula log10(\code{min.ratio}*\code{lambdamax}).
##'	}
##' }
##'
##' @return An object of class "cv.fa" for which a \code{plot} method
##' is available.
##'
##' @seealso \code{\linkS4class{fusedanova}}, \code{\link{plot.cv.fa}}
##' and \code{\linkS4class{cv.fa}}.
##'
##' @examples \dontrun{
##' data(aves)
##' cv.out <- cv.fa(aves$weight, aves$family)
##' V100.cv.out <- cv.fa(aves$weight, aves$family, V=50)
##' }
##'
##' @keywords models, regression
##' @name cv.fa
##' @aliases cv.fa
##' @rdname cv.fa
##'
##' @export
cv.fa <- function(x,
                  class,
                  K = 10,
                  folds= split(sample(1:length(class)), rep(1:K, length=length(class))),
                  lambdalist = NULL,
                  V = 1,
                  verbose = TRUE,
                  ...) {

	## =============================================================
	## INITIALIZATION & PARAMETERS RECOVERY

	user <- list(...)
	defs <- default.args.cv()
	args <- modifyList(defs, user)

	if (Sys.info()[['sysname']] == "Windows") {
		args$mc.cores <- 1 # Windows does not support fork
	}

        if(!inherits(x, c("matrix", "Matrix")))
          x <- as.matrix(x)
	n <- length(class)
	p <- ncol(x)

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
		if (Sys.info()[['sysname']] == "Windows") {
			if(verbose){
                          warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
			}
		}
		args$checkargs = FALSE # not to check again in fused anova first call
	}

	## =============================================================
	# FIRST RUN ON THE ALL DATASET
	#  get the lambda max and a fused anova object on which we can use predict later

	firstrun = do.call(fusedanova,c(list(x=x,class=class),args))
	lambdamax = max(unlist(lapply(firstrun@result,function(x){max(x$table[,"lambda"])})))

	# generating the grid of lambda by ascending order
	if (is.null(lambdalist)) {
		if (args$log.scale == FALSE){
			args$lambdalist = seq(0,lambdamax,len=args$nlambda)
		}else{
			args$lambdalist = 10 ^ seq(log10(args$min.ratio*lambdamax),log10(lambdamax), len=args$nlambda)
		}
	}else{
		args$lambdalist = sort(lambdalist,decreasing=FALSE)
	}

	# overwrite nlambda and K
	args$nlambda = length(args$lambdalist)
	K <- length(folds)

	# list format for x and factor for y
	x = split(x, rep(1:p, each = n))
	if(!is.factor(class)){class=as.factor(class)}

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
	} else{
		algoType = "With possible Splits"
	}

	one.fold <- function(k,z) {
		omit <- folds[[k]]
		if (args$standardize==TRUE){z = normalize.cv(x=z,group=class,omit=omit)}
		err <- simplecv(xtrain=z[-omit], ytrain=class[-omit], xtest = z[omit],ytest = class[omit], args=args)
		return(err = err)
	}

        if (V > 1 & verbose) {
          cat("\nAveraging ", V," ",K,"-folds cross-validation... might take a while!", sep="")
          cat("\nV =")
        } else {
          cat("\n",K,"-folds cross-validation...", sep="")
        }

        cv.list <- list()
        global.list <- list()
        folds.list <- list()
        for (v in 1:V) {
          if (V > 1 & verbose) {
            cat(" ", v)
          }

          folds.list[[v]] <- folds
          ## Errors contains a list of tables (rows <-> folds and col<-> lambdas) for each variable
          Errors <- lapply(x,function(z){
            err  <- do.call(rbind,mclapply(1:K, function(k){one.fold(k,z)}, mc.cores= args$mc.cores,
                                           mc.preschedule=ifelse(K > 10,TRUE,FALSE)))
		return(err)
          })

                                        # CV return the mean of errors per fold and the std error.
          CV <- lapply(Errors,function(err){
            err2 = err^2
            meanerr <- colMeans(err)
            meanerr2 <- colMeans(err2)
            meanerr2 <- sqrt(1/(K-1)*(meanerr2 - meanerr^2)) # erreur sur la moyenne
            lambda.min = max(args$lambdalist[meanerr <= min(meanerr)])
                                        # if several lambda.min, we take the higher lambda
            return(list(cv.error = data.frame(err=meanerr,sd=meanerr2,lambdalist = args$lambdalist),
                        lambda.min=lambda.min))
          })

          if (p>1){
            err = aaply(laply(Errors,as.matrix),c(2,3),mean) # mean on all variables for each fold
            err2 = err^2
            meanerr <- colMeans(err)
            meanerr2 <- colMeans(err2)
            meanerr2 <- sqrt(1/(K-1)*(meanerr2 - meanerr^2))
            lambda.min = max(args$lambdalist[meanerr <= min(meanerr)])
            global = list(cv.error = data.frame(err=meanerr,sd=meanerr2,
                            lambdalist = args$lambdalist),lambda.min=lambda.min)
          }else{
            global=CV[[1]]
          }
          global.list[[v]] <- global
          cv.list[[v]] <- CV

          ## new folds
          folds <- split(sample(1:length(class)), rep(1:K, length=length(class)))
        }

        cv <- as.data.frame(apply(sapply(global.list, function(l) { l$cv.error}), 1, simplify2array))

        lb <- mean(sapply(global.list, function(l) { l$lambda.min}))

        cv <- data.frame(err=tapply(cv$err, cv$lambdalist, mean),
                         sd =tapply(cv$sd, cv$lambdalist, function(x) sqrt(mean(x^2)/V)),
                         lambdalist = args$lambdalist)
        rownames(cv) <- 1:nrow(cv)
        global <- list(cv.error= cv, lambda.min=lb)

        CV <- lapply(1:p, function(i) {
          cv <- as.data.frame(apply(sapply(cv.list, function(l) { l[[i]]$cv.error}), 1, simplify2array))
          lb <- mean(sapply(cv.list, function(l) { l[[i]]$lambda.min}))
          cv <- data.frame(err=tapply(cv$err, cv$lambdalist, mean),
                           sd =tapply(cv$sd, cv$lambdalist, function(x) sqrt(mean(x^2)/V)),
                           lambdalist = args$lambdalist)
          rownames(cv) <- 1:nrow(cv)
          return(list(cv.error= cv, lambda.min=lb))
        })

	return(new("cv.fa",
			byvariable = CV,
			global = global,
			folds  = ifelse(V == 1, folds, folds.list),
			lambdalist =args$lambdalist,
			algorithm = algoType))

}

# return error
simplecv<-function(xtrain,ytrain,xtest,ytest,args){

	# Objective : xm and xmtest should have the same length for c++ code
	# we modify ytest as a factor containing the levels of y and use tapply
	# if a level is not present in ytest, we replace NA by 0 so that it works and the error calculation is still correct
	# if a level is in ytest but not in y it's discarded
	ytrain = factor(ytrain) # (drop useless factor)

	index = ytest %in% ytrain # what to discard.
	ytest = factor(ytest[index],levels=levels(ytrain)) # discard but keep the ytrain levels
	xtest =xtest[index]  # same

	ngroup = tapply(ytrain,ytrain,length) # vector of number by group
	xm = tapply(xtrain,ytrain,mean)
	xv = rep(0,length(xm))
	if (args$weights %in% varNeededWeights){ # var needed if weights are of welch or ttest type
		xv = tapply(xtrain,ytrain,var)
		xv[is.na(xv)] <- 0
	}

	ngrouptest = tapply(ytest,ytest,length) # vector of number by group
	ngrouptest[is.na(ngrouptest)]=0

	xmtest = tapply(xtest, ytest, mean)
	xmtest[is.na(xmtest)] <- 0

	xvtest = tapply(xtest,ytest,var)*(ngrouptest-1)
	xvtest[is.na(xvtest)] <- 0

	# sort from the smallest beta to the highest
	ordre = order(xm)
	xm = xm[ordre]
	ngroup =ngroup[ordre]
	ngrouptest = ngrouptest[ordre]
	xv =xv[ordre]
	xmtest = xmtest[ordre] # sous forme de data frame plus rapide ????

	# decomposition of error (Huygens):
	# \sum_i{(\hat{Y_i}(\lambda)-Y_i)^2} = sum_k{ngroup(k)*(\hat{Y}_k - sum(Y_i in k))} + sum{Var(group_k)}
	errVar = sum(xvtest)

	if (args$splits==1){
          errEst  <- .Call("noSplitcv",R_x=xm,R_xv=xv, R_ngroup=ngroup, R_xtest =xmtest, R_ngrouptest=ngrouptest, R_args=args, PACKAGE="fusedanova")
	}else{
          errEst  <- .Call("Splitcv",R_x=xm,R_xv=xv, R_ngroup=ngroup, R_xtest =xmtest, R_ngrouptest=ngrouptest, R_args=args, PACKAGE="fusedanova")
	}

	errEst = errEst + errVar

	return(errEst)
}

###################################
# normalization of a vector for cv
###################################
normalize.cv <- function(x,group,omit){
	xtrain = x[-omit]
	grouptrain = group[-omit]
	Pooled = (nlevels(grouptrain)!= length(xtrain))
	if (Pooled ==FALSE){
		res = (x-mean(xtrain))/as.numeric((sqrt(var(xtrain))))
	}else{
		ntrain = length(xtrain)
		m  = mean(xtrain)
                ngrouptrain = tabulate(grouptrain)
                ngrouptrain = ngrouptrain[ngrouptrain != 0]
		s = rowsum(xtrain^2,grouptrain) - (1/ngrouptrain)*(rowsum(xtrain,grouptrain))^2
		s[which(ngrouptrain==1)]=0
		if (sum(s)==0){
			res = (x-mean(xtrain))/as.numeric((sqrt(var(xtrain))))
		}else{
			s = sqrt(sum(s)/(ntrain-length(ngrouptrain)))
			res = (x-m)/s
		}
	}
	return(res)
}

#############################
# default args
#############################
default.args.cv <- function() {
	return(list(
		weights = "default",
		W=matrix(nrow=0, ncol=0),
		gamma = 0 ,
		standardize = TRUE,
		splits = 0,
		epsilon =10^-10,
		checkargs = TRUE,
		nlambda =100,
		log.scale = TRUE,
		min.ratio = 1e-8,
		mc.cores = detectCores(),
		verbose = FALSE,
		mxSplitSize = 100
    ))
}

# constant list of treated weights
possibleWeights = c("default","laplace","gaussian","adaptive","naivettest","ttest","welch","personal")
varNeededWeights = c("welch", "naivettest", "ttest")
