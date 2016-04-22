##' Fit a Fused ANOVA model
##'
##' Adjust a penalized ANOVA model with Fused-LASSO (or Total
##' Variation) penality, ie. a sum of weighted \eqn{\ell_1}{l1}-norm
##' on the difference of each coefficient. See details below.
##'
##' @param x matrix (or a numeric vector) which rows represent
##' individuals and columns independant variables.
##'
##' @param class vector or factor giving the class of each
##' individual. If missing, \code{1:nrow(x)} is used (clustering mode
##' with one individual per class).
##'
##' @param ... list of additional parameters to overwrite the defaults of the
##' fitting procedure. Include :
##' \itemize{%
##'
##' \item{\code{weights}: } {character; which type of weights is
##' supposed to be used.  The supported weights are :
##' \code{"default"}, \code{"laplace"}, \code{"gaussian"},
##' \code{"adaptive"}, \code{"naivettest"}, \code{"ttest"},
##' \code{"welch"} and \code{"personal"}. Se details below.  By
##' default, its value is \code{"default"}.}
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
##' The optimization problem solved by fused-ANOVA is
##' \if{latex}{\deqn{%
##' \hat{\beta}_{\lambda} = \arg \min_{\beta}
##' \left\{\sum_{k=1}^K \sum_{i=1}^{n_k} \left(Y_{ik}-\beta_k \right)^2
##' + \lambda \sum_{k,\ell} w_{kl} \left|\beta_k - \beta_\ell \right|\right\}}}
##' \if{html}{\out{ <center> &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub></sub> =
##' argmin<sub>&beta;</sub> sum<sub>k</sub> sum_i (Y<sub>ik</sub> - &beta<sub>k</sub>)<sup>2</sup>
##' + &lambda; sum<sub>k,l</sub> w<sub>k,l</sub>
##' &#124; &beta;<sub>k</sub> - &beta;<sub>l</sub> &#124;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda) = argmin_beta sum_k sum_i (Y_ik - beta_k)^2
##' + lambda sum_k sum_l w_kl | beta_k - beta_l|,}}
##'
##' where \eqn{Y_{ik}}{Y_ik} is the intensity of a continuous random
##' variable for sample \eqn{i}{i} in condition \eqn{k}{k} and
##' \eqn{\beta_k}{beta_k} is the mean parameter of condition
##' \eqn{k}{k}. We denote by \eqn{K}{K} the total number of conditions
##' and \eqn{n_k}{n_k} the number of sample in each condition.
##'
##' More details related to the weights are coming...
##'
##' @seealso See also \code{\linkS4class{fusedanova}},
##' \code{\link{plot,fusedanova-method}} and \code{\link{cv.fa}}.
##' @name fusedanova
##' @rdname fusedanova
##' @keywords models, regression
##'
##' @examples \dontrun{
##' data(aves)
##' fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
##' plot(fa.laplace, labels=aves$order)
##'
##' fa.ttest <- fusedanova(x=aves$weight, class=aves$family, weights="naivettest")
##' plot(fa.ttest, labels=aves$order)
##'
##' fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
##' plot(fa.ada, labels=aves$order)
##' }
##'
##' @export
fusedanova <- function(x, class,...) {
    
    ## ===================================================
    ## GETTING DEFAULT PARAMS
    user <- list(...)
    defs <- default.args.fa()
    args <- modifyList(defs, user)
    
    ## ====================================================
    ## BASIC CONVERSIONS
    if (is.matrix(x)) {
        p <- ncol(x); n <- nrow(x) 
    } else {
        p <- 1; n <- length(x)
    }
    if (missing(class)) {
        class <- factor(1:n)
    } else if (!is.factor(class)){
        class <- as.factor(class)
    }
    
    ## ===================================================
    ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
    if (args$checkargs) {
        if (any(is.na(x)))
            stop("NA value in x not allowed.")
        if (n != length(class))
            stop("x and y have not correct dimensions")
        if (nlevels(class)==1)
            stop("y has only one level.")
        if(!(args$weights %in% possibleWeights))
            stop("Unknown weight parameter formulation. Aborting.")
        if (Sys.info()[['sysname']] == "Windows") {
            if(args$verbose){
                warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
            }
            args$mc.cores <- 1
        }
        if (length(args$lambdalist)!=0){ # sort by descending for effective prediction
            args$lambdalist = sort(args$lambdalist, decreasing = TRUE)
        }    
    }
    
    ## ====================================================
    ## SPLIT OCCURENCE TESTING AND LAUNCH
    if (args$splits==0){ # if we let the programm choose, overwrite the value
        if (args$weights=="default"||(args$weights %in% c("laplace","gaussian","adaptive") && args$gamma>=0) ||length(unique(class))<3){
            args$splits <- 1
        }else{
            args$split <- 2
        }
    }
    if (args$split==1){
        algoType <- "No Split"
        if(args$verbose){cat("\nPath calculated without split")}
    } else {
        algoType <- "With possible Splits"
        if(args$verbose){cat("\nPath calculated with possible splits")}
    }
    
    if (p > 1) {
        ## multivariate case (distributed)
        x <- split(x, rep(1:p, each = n)) 
        if (args$standardize == TRUE) # normalization
            x <- lapply(x,function(z){normalize(z,class)})
        result <- mclapply(x,function(z) {calculatepath(z,class,args)},
                           mc.cores= args$mc.cores, mc.preschedule=ifelse(p > 100, TRUE, FALSE)) # path algorithm
    } else {
        ## univariate case
        if (args$standardize==TRUE) # normalization
            x <- normalize(x,class)
        result <- list(calculatepath(x,class,args))
    }

    ## formating the output
    if (length(args$lambdalist) == 0){ # default parameter
        l <- numeric(0)
        prediction =  list()
    } else {
        l <- args$lambdalist
        prediction <- lapply(result, function(z){list(table=z$prediction, order = z$order)})
    }    
    result <- lapply(result, function(z){list(table=z$table, order = z$order)})

    ## small warning on last beta
    if (args$standardize==TRUE && abs(result[[1]]$table[1,1])>10^(-8)){
        warning("There may be some approximation errors (Beta(lambda_max) far from 0). You may want to lower the gamma if you are using one.")}
    
    return(new("fusedanova",result=result, prediction=prediction, classes = class, weights=args$weights, lambdalist=l,algorithm=algoType))
}


#########################################
# Calculate path with or without splits #
#########################################
# calculate the path for one gene
# x the vector of data, group the vector group belonging, args all the rest
calculatepath <- function(x, group, args){
  
    ngroup <- tabulate(group)        # vector of number by group
    xm     <- rowsum(x,group)/ngroup # vector of empirical means
    ordre  <- order(xm)              # vector of empirical means
    xm     <- xm[ordre]              # sorted version
    ngroup <- ngroup[ordre]          # same for the group sizes

    ## this additional variable is needed for "Welch" or "ttest" weights
    xv <- numeric(0) 
    if (args$weights %in% varNeededWeights){        
        xv <- ngroup/(ngroup-1)*(rowsum(x^2,group)[ordre]/ngroup - xm^2)
    }
        
    if (args$splits==1){
        res  <- .Call("noSplit",R_x=xm,R_xv=xv,R_ngroup=ngroup, R_args=args, PACKAGE="fusedanova2")
    } else{
        res  <- .Call("withSplit",R_x=xm,R_xv=xv,R_ngroup=ngroup, R_args=args, PACKAGE="fusedanova2")
    }
    
    return(list(table=res$res, prediction=res$pred, order=ordre))
  
}

#############################
# normalization of a vector
#############################
normalize <- function(x,group){
  Pooled = (nlevels(group)!= length(x))
  if (Pooled == FALSE){
    res = (x-mean(x))/as.numeric((sqrt(var(x))))
  } else {
    n <- length(x)
    m <- mean(x)
    ngroup <- tabulate(group)
    s <- 1/(ngroup-1)*(rowsum(x^2,group) - (1/ngroup)*(rowsum(x,group))^2)
    s[which(ngroup==1)] <- 0
    if (sum(s)==0){
      res <- (x-mean(x))/as.numeric((sqrt(var(x))))
    } else{
      s <- sqrt(sum(s*(ngroup-1))/(n-length(ngroup)))
      res <- (x-m)/s
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

## constant list of treated weights
possibleWeights = c("default","laplace","gaussian","adaptive","naivettest","ttest","welch","personal")
varNeededWeights = c("welch", "ttest")
