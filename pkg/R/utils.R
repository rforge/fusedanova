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

# convert a res for 1 gene to the matrix of betas, rows : class, col : lambda
# nb supposed that we have a value for each tuple (class,beta) (only the case when prediction activated)
toMat <- function(x, class= NULL){
	nlambda = length(unique(x[,"lambda"]))
	x = x[with(x, order(x[,"lambda"], x[,"class"])), "beta"]
	x = matrix(x, ncol=nlambda, byrow = FALSE) 
	if (! is.null(class)){
		x = x[class,] 
	}
	return(x)
}

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
  require(grid)
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


