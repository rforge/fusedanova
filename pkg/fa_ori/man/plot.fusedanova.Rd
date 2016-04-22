\docType{methods}
\name{plot,fusedanova-method}
\alias{plot,fusedanova-method}
\title{Plot method for a fusedanova object}
\usage{
\\S4method{plot}{fusedanova}(x, y,
main=paste("Regularization path for variable",varvect),
xlab = expression(paste("location in the regularization path  ",lambda)),
ylab = expression(paste("optimal coefficient  ",beta)),
log.scale = TRUE, reverse   = FALSE, labels    = NULL ,
varvect = sample(1:length(slot(x,"result")),1), plot = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{fusedanova}.}

  \item{y}{used for S4 compatibility.}

  \item{main}{the main title, with a hopefully appropriate
  default definition.}

  \item{xlab}{character or expression (or a "grob") giving
  label(s) for the x-axis.}

  \item{ylab}{character or expression (or "grob") giving
  label for the y-axis.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used. Default is \code{TRUE}.}

  \item{reverse}{logical; should the X-axis be reversed?
  Default is \code{FALSE}.}

  \item{labels}{a vector of factor with labels associated
  to n samples classed by fused-ANOVA. Default is
  \code{NULL}.}

  \item{varvect}{a vector with the variables whose path
  will be plot. Default picks one at random.}

  \item{plot}{logical; indicates if the graph should be
  plotted on call. Default is \code{TRUE}.}

  \item{...}{used for S4 compatibility.}
}
\value{
a \pkg{ggplot2} object which can be plotted via the
\code{print} method.
}
\description{
Produce a plot of the solution path of a \code{fusedanova}
fit.
}
\examples{
data(aves)
fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
plot(fa.laplace, labels=aves$order)
}
\seealso{
\code{\linkS4class{fusedanova}}.

\code{\linkS4class{fusedanova}}.

\code{\linkS4class{fusedanova}}.
}

