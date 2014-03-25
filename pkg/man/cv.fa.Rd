\name{cv.fa}
\alias{cv.fa}
\title{Cross-validation function for fusedanova method.}
\usage{
cv.fa(x, class, K = 10, folds = split(sample(1:length(class)), rep(1:K,
  length = length(class))), lambdalist = NULL, V = 1, verbose = TRUE, ...)
}
\arguments{
  \item{x}{matrix whose rows represent individuals and
  columns independant variables.}

  \item{class}{vector or factor giving the class of each
  individual.}

  \item{K}{integer indicating the number of folds. Default
  is 10.}

  \item{V}{integer indicating the number of times the K
  folds CV will be averaged. Default is 1 (no averaging).}

  \item{verbose}{boolean for verbose mode. Default is
  \code{TRUE}.}

  \item{folds}{list of \code{K} vectors that describes the
  folds to use for the cross-validation. By default, the
  folds are randomly sampled with the specified K. The same
  folds are used for each values of \code{lambdalist}.}

  \item{lambdalist}{list of \code{lambda} penalty
  parameters used in the cross validation. By default,
  lambdalist is NULL and calculated using the maximum
  \code{lambda} and the parameter \code{nlambda}.}

  \item{...}{list of additional parameters to overwrite the
  defaults of the fitting procedure. See the corresponding
  documentation (\code{\link{fusedanova}}). Also include :
  \itemize{%

  \item{\code{nlambda}: } {integer; the length of the
  \code{lambdalist} vector, by default 100}

  \item{\code{log.scale}: } {boolean; should a logarithmic
  scale be used during the creation of \code{lambdalist}.
  By default, FALSE. }

  \item{\code{min.ratio} : }{ numeric parameter setting the
  smallest value of lambdalist in the \code{log.scale} case
  with the formula
  log10(\code{min.ratio}*\code{lambdamax}).  } }}
}
\value{
An object of class "cv.fa" for which a \code{plot} method
is available.
}
\description{
Function that computes K-fold cross-validated error of a
\code{fusedanova} fit. Possibility to perform V times the
whole K-fold CV to average the error and reduce the reduce
the CV standard error around the minimum.
}
\examples{
\dontrun{
data(aves)
cv.out <- cv.fa(aves$weight, aves$family)
V100.cv.out <- cv.fa(aves$weight, aves$family, V=50)
}
}
\seealso{
\code{\linkS4class{fusedanova}},
\code{\link{plot,cv.fa-method}} and
\code{\linkS4class{cv.fa}}.
}
\keyword{models,}
\keyword{regression}

