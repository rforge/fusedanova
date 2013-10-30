\name{fusedanova}
\alias{fusedanova}
\title{Fit a Fused ANOVA model}
\usage{
  fusedanova(x, class, ...)
}
\arguments{
  \item{x}{matrix which rows represent individuals and
  columns independant variables.}

  \item{class}{vector or factor giving the class of each
  individual.}

  \item{...}{list of additional parameters to overwrite the
  defaults of the fitting procedure. Include : \itemize{%

  \item{\code{weights}: } {character; which type of weights
  is supposed to be used. The supported weights are :
  \code{"default"}, \code{"laplace"}, \code{"gaussian"},
  \code{"adaptive"}, \code{"naivettest"}, \code{"ttest"},
  \code{"welch"} and \code{"personal"}. By default, its
  value is \code{"default"}.}

  \item{\code{gamma}: } {numeric; the \eqn{\gamma}{gamma}
  parameter needed for \code{"laplace"}, \code{"gaussian"}
  and \code{"adaptive"} weights. By default, 0.}

  \item{\code{W}: } {numeric matrix; the matrix of weights
  needed if the \code{"personal"} weights were selected. By
  default, a matrix with zero row and zero column.}

  \item{\code{standardize}: }{ logical; should each
  variable be standardized before the calculus ? By
  default, \code{TRUE}. }

  \item{\code{splits}: }{ integer; coding for forcing split
  or nosplit algorithms : \itemize{%

  \item{\code{0} : }{Default value, let the programm decide
  which algorithm to use depending on the choosen
  \code{weights}.}

  \item{\code{1} : }{Forces the algorithm not to take the
  splits into account.}

  \item{\code{2} : }{Forces the algorithm to take the
  splits into account even if the solution paths contains
  no split.} } Note : For the moment, only the no split
  algorithm has been coded. Please ensure that your weights
  choice leads to the no split algorithm or set
  \code{split} to \code{1}. }

  \item{\code{epsilon}: }{numeric; tolerance parameter for
  numeric calculations. By default, \eqn{10^-10}{eps}. Note
  : this is currently not used.}

  \item{\code{checkargs}: }{logical; should arguments be
  checked to (hopefully) avoid internal crashes? Default is
  \code{TRUE}. Automatically set to \code{FALSE} when a
  call is made from cross-validation}

  \item{\code{lambdalist}: }{numeric vector; a set of
  \eqn{\lambda}{lambda} value for which a prediction is
  asked. By default, a null vector. If the length of
  \code{lambdalist} is not \code{0}, the \code{fusedanova}
  class returnded by the \code{fusedanova} function will
  have a not null attribute \code{prediction}. }

  \item{\code{mc.cores}: } {integer; the number of cores to
  use. The default uses all the cores available. }

  \item{\code{verbose}: } {boolean; should the code print
  out its progress. By default, FALSE. }

  \item{\code{mxSplitSize}: } {integer; the maximum size
  for a group for being checked the cores available. }

  }}
}
\value{
  an object with class \code{fusedanova}, see the
  documentation page \code{\linkS4class{fusedanova}} for
  details.
}
\description{
  Adjust a penalized ANOVA model with Fused-LASSO (or Total
  Variation) penality, ie. a sum of weighted
  \eqn{\ell_1}{l1}-norm on the difference of each
  coefficient.
}
\note{
  The optimized criterion is the following:
}
\examples{
data(aves)
fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
plot(fa.laplace, labels=aves$order)

fa.ttest <- fusedanova(x=aves$weight, class=aves$family, weights="naivettest")
plot(fa.ttest, labels=aves$order)

fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
plot(fa.ada, labels=aves$order)
}
\seealso{
  See also \code{\linkS4class{fusedanova}},
  \code{\link{plot.fusedanova}} and \code{\link{crossval}}.
}
\keyword{models,}
\keyword{regression}

