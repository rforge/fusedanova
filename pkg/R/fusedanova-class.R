##' Class "fusedanova"
##'
##' Class of object returned by the fusedanova function.
##'
##' @section Slots: \describe{
##'
##' \item{\code{result}: }{a list of objects, one per variable. Objects are lists containing :
##' \itemize{%
##' \item{\code{table}: }{a dataframe containing the \code{lambda}; \code{Beta}, value of coefficient
##' ; \code{slope}; \code{idown} the minimal index of class in the group and \code{iup}, the maximum one;
##'  for each fusion time \code{lambda}.}
##'
##' \item{\code{order}: }{a vector giving the order of means before any fuse. 
##' Needed with \code{classes} for data conversion to a more understandable format.}
##'
##' } 
##' }
##'
##' \item{\code{classes}: }{the intitulate of the classes as entered in the \code{class}
##' parameter in the \code{fusedanova} function.}
##'
##' \item{\code{prediction}: }{if a \code{lambdalist} was given to the \code{fusedanova} function,
##' a list contaning for each variable an object with the same format as the ones in 
##' \code{result}. Else, an empty list.}
##'
##' \item{\code{lambdalist}: }{the \code{lambdalist} vector parameter given in the
##'  \code{fusedanova} function. If not given, empty vector.}
##'
##' } 
##'
##' @section Methods:
##' Specific plotting, predict and conversion methods are available and documented
##' (\code{\link{plot.fusedanova}}, \code{\link{predict.fusedanova}}, 
##' \code{\link{dataconvert.fusedanova}}).
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot.fusedanova}}, \code{\link{predict.fusedanova}}
##' \code{\link{dataconvert.fusedanova}} and \code{\link{fusedanova}}.
##'
##' @name fusedanova-class
##' @rdname fusedanova-class
##'
##' @exportClass fusedanova
##' @exportMethod predict
##' @exportMethod dataconvert
setClass("fusedanova",
	representation = representation(
		result  = "list", # list of triple (dataframe,vector,vector)
		classes = "factor",
		prediction = "list",
		lambdalist = "vector",
		algorithm = "character"
    )
)

##' Yeah well I need a generic, what do I comment ?
setGeneric ( name= "dataconvert",
	def = function (object,...){ standardGeneric ("dataconvert")}
)

##' Conversion method for a fusedanova object
##'
##' Produce for each object (dataframe, vector of oreder, vector of class) in the list of 
##' \code{result} or \code{prediction} (if \code{predicted} == \code{TRUE}) a new dataframe
##' containing for each fusion the \code{lambda}, \code{beta}, \code{slope}, \code{class}.
##'
##' @param object an object of class \code{fusedanova}.
##' @param predicted logical; if \code{TRUE}, the return is calculated on the \code{result} slot,
##'  else on the \code{prediction} slot. By default, \code{FALSE}.
##' @param listformat logical; does the return should be a list of dataframes for each variable 
##' or only one dataframe with one more column indexing the variables. By default, \code{FALSE}.
##' @param ... used for S4 compatibility.
##' 
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @name dataconvert,fusedanova-method
##' @aliases dataconvert,fusedanova-method
##' @aliases dataconvert.fusedanova
##' @docType methods
##' @rdname dataconvert.fusedanova
##'
##' @examples \dontrun{
##' ## Not done yet
##' }
##'
##' @export
setMethod("dataconvert", "fusedanova", function (object, predicted=FALSE, formating = c("df","list")[1], labels =FALSE,...){
	# set the good format for plot or visualization
	# for the moment extremely dirty
	new = NULL
	if (predicted==FALSE){
		res=object@result
	}else {
		res = object@prediction
	}
    
	class = object@classes
	
	if (formating == "df"){ # dataframe with vars column give the variable
		if (object@algorithm == "No Split"){
			for (l in 1:length(res)){ 
				x  = res[[l]]
				ordre = x$order
				x =x$table
				for (i in 1:nrow(x)){
					if (x$iup[i]>=x$idown[i]){
						for (j in 1:(x$iup[i]-x$idown[i]+1)){
							new = rbind(new, c(x$beta[i],x$lambda[i],x$slope[i],ordre[(x$idown[i]:x$iup[i])[j]],l))
						}
					} 
				}
			}

		}else{
			for (l in 1:length(res)){ 
				x  = res[[l]]
				ordre = x$order
				x =x$table
				x$class = ordre[x$class]
				new = rbind(new,cbind(x,rep(l,nrow(x))))
			} 
		}
			
		new = as.data.frame(new)
		colnames(new) =c("beta","lambda","slope","class","vars")
		if (labels){new$class = as.factor(levels(class)[new$class])}
		new =unique(new[order(new[,"vars"], -new[,"lambda"], new[,"class"]), ])

	}else if(formating == "list"){ # list : 1 element of list = one variable
	
		if (object@algorithm == "No Split"){
			new = lapply(res, function(z){
				new = NULL
				ordre = z$order
				x =z$table
				for (i in 1:nrow(x)){
					if (x$iup[i]>=x$idown[i]){
						for (j in 1:(x$iup[i]-x$idown[i]+1)){
							new = rbind(new, c(x$beta[i],x$lambda[i],x$slope[i],ordre[(x$iup[i]:x$idown[i])[j]]))
						}
					} 
				}
				new = as.data.frame(new)
				colnames(new) =c("beta","lambda","slope","class")
				if (labels){new$class = as.factor(levels(class)[new$class])}
				new = unique(new[order(-new[,"lambda"], new[,"class"]), ])		
			})
			
		}else{
			new = lapply(res, function(z){
					new = NULL
					ordre = z$order
					new =z$table
					new$class = ordre[x$class]
					if (labels){new$class = as.factor(levels(class)[new$class])}
					new = unique(new[order(-new[,"lambda"], new[,"class"]), ])
			})
		}	
	} 
	
	# matrix format, cost more, to do ???. <-----------------------

	return(new)
})


#setMethod("print", "fusedanova", definition =
#   function(object) {print(object)}
#)

#setMethod("show", "fusedanova", definition =
#   function(object) {print(object)}
#)


##' Plot method for a fusedanova object
##'
##' Produce a plot of the solution path of a \code{fusedanova} fit.
##' 
##' @param x an object of class \code{fusedanova}. 
##' @param y used for S4 compatibility.
##' @param common type of plot. By default \code{"l"}.
##' @param main the main title, with a hopefully appropriate default definition.
##' @param xlab character or expression (or a "grob") giving label(s) for the x-axis.
##' @param ylab character or expression (or "grob") giving label for the y-axis.
##' @param strip logical flag or function. If FALSE, strips are not drawn.
##' @param layout the layout of the plot, \code{NULL} by default.
##' @param ... used for S4 compatibility.
##' 
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @name plot,fusedanova-method
##' @aliases plot,fusedanova-method
##' @aliases plot.fusedanova
##' @docType methods
##' @rdname plot.fusedanova
##'
##' @examples \dontrun{
##' ## Not done yet
##' }
##'
##' @export
setMethod("plot", "fusedanova", definition =
	function(x,y,type="l",
	main="The entire regularization path of optimal solutions for each variable", # main title
	xlab=expression(paste("location in the regularization path  ",lambda)), # horizontal axis
	ylab=expression(paste("optimal coefficient  ",beta)), #vertical axis
	log.scale = FALSE,
 	strip=strip.custom(strip.names=TRUE),
	layout = NULL, #layout for the plot
	varvect = NULL,
	...){

	df  = dataconvert(x)

	if (log.scale) {
		epsilon = 10^-1
		lambdalist = unique(df$lambda[df$lambda!=0])
		lambdalist = c(min(lambdalist)*epsilon, lambdalist)
		df = predict(x,lambda=lambdalist)
		df$lambda <- log10(df$lambda)
	}else{
	
	}
	
	if (!is.null(varvect)){
		df <- unique(subset(df, vars %in% varvect)[,c("beta","lambda","slope","class","vars")])
	}

	if (is.null(layout)){
		layout=c(1,nlevels(as.factor(df$vars)))
	}
	
	xyplot(beta~lambda|vars,df,group=class,type=type,layout=layout,
			xlab=xlab,ylab=ylab,main=main,strip=strip,...)

})


##' Predict method for a fusedanova object
##'
##' Produce a prediction for a vector of \code{lambda} parameter and an array of \code{class}.
##' 
##' @param object an object of class \code{fusedanova}. 
##' @param y a vector of \code{class}. By default, \code{NULL}. If \code{NULL}, all classes
##' are predicted.
##' @param lambda a numeric vector giving the list of \eqn{\lambda}{lambda} for which to predict.
##' By default, \code{NULL}. If \code{NULL}, it is set to the \code{lambdalist} slot 
##' of \code{object}. If this slot is empty, \code{lambda} is set to the fusion times detected in
##' the \code{fusedanova} function.
##'
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @name predict,fusedanova-method
##' @aliases predict,fusedanova-method
##' @aliases predict.fusedanova
##' @docType methods
##' @rdname predict.fusedanova
##'
##' @examples \dontrun{
##' ## Not done yet
##' }
##'
##' @export
setMethod("predict", "fusedanova", definition =
	function (object, y= NULL, lambda=NULL, labels = FALSE)  {
	require(plyr) 	
	if (is.null(lambda)){ # no new grid
		if (length(object@lambdalist)==0){ # no pred was asked when launching fused anova
			d=dataconvert(object,labels =labels)
			if (!is.null(y)){
				d=subset(d,class %in% y)
			}
		}else{
			d=dataconvert(object,predicted=TRUE,labels =labels)
			if (!is.null(y)){
				d=subset(d,class %in% y)
			}
		}
		
	}else{ # linear interpolation 
	  
		res = dataconvert(object,labels=labels)
		calc1 <- function(d){
			dm <- merge(d,data.frame(lambda),all=TRUE)
			tofill <- transform(dm,class=class[1],vars=vars[1],slope=slope[1])
			fillin(tofill,"beta")
		}

		d <- ddply(res,.(class,vars),calc1)
		lambdalist <- lambda
		d <- unique(subset(d,lambda %in% lambdalist)[,c("beta","lambda","slope","class","vars")])
	   
		if (!is.null(y)){ # no new data
			d=subset(d,class %in% y)
		} 
	}
	d = d[order(-d[,"lambda"], d[,"class"]),]
	return(d)
})

