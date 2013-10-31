##' Class "fusedanova"
##'
##' Class of object returned by the \code{fusedanova} function.
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
##' \item{\code{weights}: }{the weights used t eprform the fit.}
##'
##' ' \item{\code{classes}: }{the intitulate of the classes as entered in the \code{class}
##' parameter in the \code{fusedanova} function.}
##'
##' \item{\code{prediction}: }{if a \code{lambdalist} was given to the \code{fusedanova} function,
##' a list contaning for each variable an object with the same format as the ones in
##' \code{result}. Else, an empty list.}
##'
##' \item{\code{lambdalist}: }{the \code{lambdalist} vector parameter given in the
##'  \code{fusedanova} function. If not given, empty vector.}
##'
##' \item{\code{algorithm}: }{a character string indicating wether splits where allowed or not.}
##'
##' }
##'
##' @section Methods:
##' Specific plotting, predict and conversion methods are available and documented
##' (\code{\link{plot.fusedanova}}, \code{\link{predict.fusedanova}},
##' \code{\link{dataconvert.fusedanova}}).
##'
##' @aliases predict,fusedanovaquadrupen-method
##' print,fusedanova-method show,fusedanova-method
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
##' @exportMethod print
##' @exportMethod show
##'
##' @importFrom stats predict
##'
setClass("fusedanova",
	representation = representation(
		result  = "list", # list of triple (dataframe,vector,vector)
		classes = "factor",
		prediction = "list",
		lambdalist = "vector",
		algorithm = "character",
                weights    = "character"
    )
)

setMethod("print", "fusedanova", definition =
   function(x, ...) {
     cat("Fused-ANOVA fit with ", x@weights, "weigths.", x@algorithm, "can occur along the path.\n")
     cat("- number of variables:", length(x@result),"\n")
     cat("- number of classes  :", nlevels(x@classes), "\n")
     invisible(x)
   }
)

setMethod("show", "fusedanova", definition =
   function(object) {print(object)}
)

##' Conversion method for a fusedanova object
##'
##' Convert the compressed fused-ANOVA path to a more handy but more
##' memory demanding data.frame format. This is typically used for
##' plotting purposes.
##'
##' @param object an object of class \code{fusedanova}.
##' @param predicted logical; if \code{TRUE}, the return value expands
##' the \code{result} slot of the original \code{fusedanova}
##' object. Otherwise, the \code{prediction} slot is
##' explanded. Default is \code{FALSE}.
##' @param listformat logical; does the return value should be a list
##' of dataframes for each variable or only one dataframe with an
##' additional column indexing the variables. By default,
##' \code{FALSE}.
##' @param ... used for S4 compatibility.
##'
##' @name dataconvert,fusedanova-method
##' @aliases dataconvert,fusedanova-method
##' @aliases dataconvert.fusedanova
##' @aliases dataconvert
##' @docType methods
##' @rdname dataconvert.fusedanova
##'
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @examples \dontrun{
##' data(aves)
##' fa <- fusedanova(x=aves$weight, class=aves$family)
##' dataconvert(fa)
##' }
##'
##' @exportMethod dataconvert
setGeneric ( name= "dataconvert",
	def = function (object,...){ standardGeneric ("dataconvert")}
)
setMethod("dataconvert", "fusedanova",
   function (object, predicted=FALSE, formating = c("df","list")[1], labels =FALSE,...){
     ## set the good format for plot or visualization
     ## for the moment extremely dirty
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
       } else {
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

     } else if(formating == "list"){ # list : 1 element of list = one variable

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

     ## matrix format, cost more, to do ???. <-----------------------

     return(new)
})


##' Plot method for a fusedanova object
##'
##' Produce a plot of the solution path of a \code{fusedanova} fit.
##'
##' @param x an object of class \code{fusedanova}.
##' @param y used for S4 compatibility.
##' @param main the main title, with a hopefully appropriate default definition.
##' @param xlab character or expression (or a "grob") giving label(s) for the x-axis.
##' @param ylab character or expression (or "grob") giving label for the y-axis.
##' @param log.scale logical; indicates if a log-scale should be used. Default is \code{TRUE}.
##' @param reverse logical; should the X-axis be reversed? Default is \code{FALSE}.
##' @param labels a vector of factor with labels associated to n
##' samples classed by fused-ANOVA. Default is \code{NULL}.
##' @param varvect a vector with the variables whose path will be plot. Default picks one at random.
##' @param ... used for S4 compatibility.
##'
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @name plot,fusedanova-method
##' @aliases plot,fusedanova-method
##' @aliases plot.fusedanova
##' @docType methods
##' @rdname plot.fusedanova
##' @seealso \code{\linkS4class{fusedanova}}.
##'
##' @examples \dontrun{
##' data(aves)
##' fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
##' plot(fa.laplace, labels=aves$order)
##' }
##'
##' @export
setMethod("plot", "fusedanova", definition =
  function(x, y,
   main  = paste("The entire regularization path of optimal solutions for variable",varvect), # main title
   xlab = expression(paste("location in the regularization path  ",lambda)), # horizontal axis
   ylab = expression(paste("optimal coefficient  ",beta)), # vertical axis
   log.scale = TRUE,
   reverse   = FALSE,
   labels    = NULL ,
   varvect   = sample(1:length(x@result),1), plot = TRUE, ...) {

  df <- dataconvert(x)

  if (log.scale) {
    epsilon = 10^-1
    lambdalist = unique(df$lambda[df$lambda!=0])
    lambdalist = c(min(lambdalist)*epsilon, lambdalist)
    df = predict(x,lambda=lambdalist)
    df$lambda <- log10(df$lambda)
  }

  if (!is.null(varvect)){
    df <- unique(subset(df, df$vars %in% varvect)[,c("beta","lambda","slope","class","vars")])
    df$vars <- as.factor(paste("var.", df$vars))
  }

  if (!is.null(labels)) {
    df$labels <- as.factor(unique(data.frame(x@classes, labels))[df$class, 2])
  } else {
    df$labels <- as.factor(df$class)
  }

  d <- ggplot(data=df,aes(x=lambda,y=beta, colour=labels, group=as.factor(class))) +
    geom_line() + labs(x=xlab, y=ylab, title=main) +
      geom_hline(yintercept=0, alpha=0.5, linetype="dotted")  +
        scale_colour_discrete(guide = guide_legend(title = "Classification"))

  if (reverse==TRUE) {
    d <- d + scale_x_reverse()
  }

  if (is.null(labels) | nlevels(df$class) > 20) {
    d <- d + theme(legend.position="none")
  }

  if (!is.null(varvect)){
    d <- d + facet_grid(.~vars)
  }

  if (plot) {print(d)}

  return(d)
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
##' @examples
##' data(aves)
##' fa <- fusedanova(x=aves$weight, class=aves$family, weight="laplace", gamma=5)
##' predict(fa, labels=aves$order)
##'
##' @export
setMethod("predict", "fusedanova", definition =
	function (object, y= NULL, lambda=NULL, labels = FALSE)  {
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

# function from Toby Hocking
# complete the alpha dataframe in the predefine col
# with linear interpolation
fillin <- function (alpha, col)
{
  y <- alpha[, col]
  na <- is.na(y)
  x <- alpha$lambda
  y[na] <- approx(x[!na], y[!na], x[na])$y
  y[is.na(y)] <- rev(y[!na])[1]
  alpha[, col] <- y
  return(alpha)
}
