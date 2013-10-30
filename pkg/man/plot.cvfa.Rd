\docType{methods}
\name{plot,cvfa-method}
\alias{plot.cvfa}
\alias{plot,cvfa-method}
\title{Plot method for cross validated error of a \code{fusedanova} model}
\arguments{
  \item{x}{output of a \code{crossval} run (must be of
  class \code{cvfa}).}

  \item{y}{used for S4 compatibility.}

  \item{varvect}{integer vector; give the index of all
  variables to plot separately. By default, varvect is NULL
  and the plot returns the plot of the \code{global}
  attribute of a \code{cvfa} instance.}

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
  a \pkg{ggplot2} object which can be plotted via the
  \code{print} method if no list of variable was given or a
  list of \pkg{ggplot2} object which can be plotted via the
  \code{multiplot} function.
}
\description{
  Produce a plot of the cross validated error of a
  \code{fusedanova} model.
}
\examples{
\dontrun{
#Not done yet
}
}

