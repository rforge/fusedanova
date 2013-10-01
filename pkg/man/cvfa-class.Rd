\docType{class}
\name{cvfa-class}
\alias{cvfa-class}
\title{Class "cvfa"}
\description{
  Class of object returned by a cross-validation performed
  through the \code{crossval} method.
}
\section{Slots}{
  \describe{ \item{\code{byvariable}:}{list (for each
  variable) of list containing a data frame (containing the
  cross validation error and its min \code{lambda}}
  \item{\code{global}:}{almost identical to byvariable
  although it contains only one list containing a dataframe
  (of the mean error across all variables and its
  associated standard error for each lambda ) and the
  lambda min} \item{\code{lambdalist}:}{vector of
  \eqn{\lambda} for which each cross-validation has been
  performed.} \item{\code{folds}:}{list of \code{K} vectors
  indicating the folds used for cross-validation.}
  \item{\code{algorithm}:}{Indicates which of the two
  algorithms was used.} }

  The specific \code{\link{plot.cvfa}} method is
  documented.
}
\seealso{
  See also \code{\link{plot.cvfa}} and
  \code{\link{crossval}}.
}
\keyword{class}

