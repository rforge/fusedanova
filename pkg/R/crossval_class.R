##' Class "cv.fa"
##'
##' Class of object returned by a cross-validation performed through
##' the \code{cv.fa} method.
##'
##' @section Slots: \describe{
##' \item{\code{byvariable}:}{list (for each variable) of list
##' containing a data frame (containing the cross validation error and its
##  associated value of standard error for each value of \code{lambda}) and the
##' min \code{lambda}}
##' \item{\code{global}:}{almost identical to byvariable although it contains
##' only one list containing a dataframe (of the mean error across all variables and
##' its associated standard error for each lambda ) and the lambda min}
##' \item{\code{lambdalist}:}{vector of \eqn{\lambda}
##' for which each cross-validation has been performed.}
##' \item{\code{folds}:}{list of \code{K} vectors indicating the folds
##' used for cross-validation.}
##' \item{\code{algorithm}:}{Indicates which of the two algorithms was used.}
##' }
##'
##' The specific \code{\link{plot}} method is documented.
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot,cv.fa-method}} and
##' \code{\link{cv.fa}}.
##'
##' @name cv.fa-class
##'
##' @rdname cv.fa-class
##'
##' @exportClass cv.fa
##'
setClass("cv.fa",
	representation = representation(
	byvariable = "list",
	global = "list",
	lambdalist = "numeric",
	folds = "list",
	algorithm = "character"
	# global and variable contains one list for each gene or one global list containing :
	# cv.error = "numeric", lambda.min = "numeric", cv.sd ="numeric" and beta.min = "numeric"
	)
)

##' Plot method for cross validated error of a \code{fusedanova} model
##'
##' Produce a plot of the cross validated error of a \code{fusedanova}
##' model.
##'
##' @usage \\S4method{plot}{cv.fa}(x, y, varvect = NULL, plotting=TRUE,
##' lambdamins=TRUE, log.scale=TRUE, reverse=FALSE,main = "Cross-validation error",...)
##' @param x output of a \code{crossval} run (must be of class
##' \code{cv.fa}).
##' @param y used for S4 compatibility.
##' @param varvect integer vector; give the index of all variables to plot separately.
##' By default, varvect is NULL and the plot returns the plot of the \code{global} attribute
##' of a \code{cv.fa} instance.
##' @param plotting logical; indicates if the graph should be plotted. Default is \code{TRUE}.
##' @param lambdamins logical; should the distribution of lambdamin be plotted ?
##' @param log.scale logical; indicates if a log-scale should be used
##' @param reverse logical; should the X-axis by reversed when \code{xvar=lambda}? Default is FALSE.  Ignored for 2D cross-validation plot.
##' @param main the main title, with a hopefully appropriate default definition.
##' @param ... used for S4 compatibility.
##'
##' @return a \pkg{ggplot2} object or a list of \pkg{ggplot2} objects
##' (if a list of variable was provided) which can be plotted via the
##' \code{print} method
##' @name plot,cv.fa-method
##' @aliases plot,cv.fa-method
##' @aliases plot.cv.fa
##' @docType methods
##' @rdname plot.cv.fa
##'
##' @seealso \code{\linkS4class{cv.fa}}.
##'
##' @examples \dontrun{
##' data(aves)
##' cv.out <- cv.fa(aves$weight, aves$family)
##' V100.cv.out <- cv.fa(aves$weight, aves$family, V=50)
##' plot(cv.out)
##' plot(V100.cv.out)
##' }
##'
##' @exportMethod plot
##' @export
setMethod("plot", "cv.fa", definition =
  function(x, y, varvect = NULL, plotting=TRUE, lambdamins=TRUE, log.scale=TRUE, reverse=FALSE,main = "Cross-validation error",...) {

    K <- length(x@folds)
    n <- length(unlist(x@folds))

	## GLOBAL CROSS-VALIDATION GRAPH
	if (is.null(varvect)){

		if (log.scale) {
			x@global$cv.error$lambdalist <- log10(x@global$cv.error$lambdalist)
		}
		d <- ggplot(x@global$cv.error, aes(x=lambdalist,y=err)) + ylab("Mean square error") + geom_point(alpha=0.3)
		d <- d + geom_smooth(aes(ymin=err-sd, ymax=err+sd), data=x@global$cv.error, alpha=0.2, stat="identity")
		if (reverse==TRUE) {
			d <- d + scale_x_reverse()
		}
		if (log.scale) {
			d <- d + xlab(expression(log[10](lambda)))
			d <- d + annotation_logticks(sides="b")
		} else {
			d <- d + xlab( expression(lambda) )
		}
		if (lambdamins & is.null(varvect)){
			if (log.scale) {
				lambda <- data.frame(xval=log10(unlist(lapply(x@byvariable, function(z){z$lambda.min}))),
                             lambda.choice=factor(1:length(x@byvariable)))
			}else{
				lambda <- data.frame(xval= unlist(lapply(x@byvariable, function(z){z$lambda.min})),
                             lambda.choice=factor(1:length(x@byvariable)))
			}
			d <- d + geom_vline(data=lambda, aes(xintercept=xval,colour=lambda.choice),
							linetype="dashed",  alpha=0.5)
			#d <- d + geom_line(data = density(lambda$xval), alpha = 0.2)
		}

		d <- d + ggtitle(main) #+

		## DO THE PLOT
		if (plotting) {print(d)}

	}else{

		# function that generate a plot of a variable
		plotunique <- function(variable,log.scale,reverse){
			if (log.scale) {
				variable$cv.error$lambdalist <- log10(variable$cv.error$lambdalist)
			}
			d <- ggplot(variable$cv.error, aes(x=lambdalist,y=err)) + ylab("Mean square error") + geom_point(alpha=0.3)
			d <- d + geom_smooth(aes(ymin=err-sd, ymax=err+sd), data=variable$cv.error, alpha=0.2, stat="identity")
			if (reverse==TRUE) {
				d <- d + scale_x_reverse()
			}
			if (log.scale) {
				d <- d + xlab(expression(log[10](lambda)))
				d <- d + annotation_logticks(sides="b")
			} else {
				d <- d + xlab(expression(lambda))
			}
			return(d)
		}

		# generate one plot for each choosen variable
		d = lapply(x@byvariable[varvect],
					function(z){plotunique(z, log.scale=log.scale,reverse=reverse)})

		if (plotting){multiplot(plotlist=d)}
	}

	return(d)
})

# Multiple plot function
# from http://www.cookbook-r.com/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(...,plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
