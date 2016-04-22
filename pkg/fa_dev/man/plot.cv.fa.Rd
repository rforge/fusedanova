\docType{methods}
\name{plot,cv.fa-method}
\alias{plot,cv.fa-method}
\title{Plot method for cross validated error of a \code{fusedanova} model}
\usage{
\\S4method{plot}{cv.fa}(x, y, varvect = NULL, plotting=TRUE,
lambdamins=TRUE, log.scale=TRUE, reverse=FALSE,main = "Cross-validation error",...)
}
\arguments{
  \item{x}{output of a \code{crossval} run (must be of
  class \code{cv.fa}).}

  \item{y}{used for S4 compatibility.}

  \item{varvect}{integer vector; give the index of all
  variables to plot separately. By default, varvect is NULL
  and the plot returns the plot of the \code{global}
  attribute of a \code{cv.fa} instance.}

  \item{plotting}{logical; indicates if the graph should be
  plotted. Default is \code{TRUE}.}

  \item{lambdamins}{logical; should the distribution of
  lambdamin be plotted ?}

  \item{log.scale}{logical; indicates if a log-scale should
  be used}

  \item{reverse}{logical; should the X-axis by reversed
  when \code{xvar=lambda}? Default is FALSE.  Ignored for
  2D cross-validation plot.}

  \item{main}{the main title, with a hopefully appropriate
  default definition.}

  \item{...}{used for S4 compatibility.}
}
\value{
a \pkg{ggplot2} object or a list of \pkg{ggplot2} objects
(if a list of variable was provided) which can be plotted
via the \code{print} method
}
\description{
Produce a plot of the cross validated error of a
\code{fusedanova} model.
}
\examples{
\dontrun{
data(aves)
cv.out <- cv.fa(aves$weight, aves$family)
V100.cv.out <- cv.fa(aves$weight, aves$family, V=50)
plot(cv.out)
plot(V100.cv.out)
}
}
\seealso{
\code{\linkS4class{cv.fa}}.
}

