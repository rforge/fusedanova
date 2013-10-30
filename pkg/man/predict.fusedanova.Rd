\docType{methods}
\name{predict,fusedanova-method}
\alias{predict.fusedanova}
\alias{predict,fusedanova-method}
\title{Predict method for a fusedanova object}
\arguments{
  \item{object}{an object of class \code{fusedanova}.}

  \item{y}{a vector of \code{class}. By default,
  \code{NULL}. If \code{NULL}, all classes are predicted.}

  \item{lambda}{a numeric vector giving the list of
  \eqn{\lambda}{lambda} for which to predict. By default,
  \code{NULL}. If \code{NULL}, it is set to the
  \code{lambdalist} slot of \code{object}. If this slot is
  empty, \code{lambda} is set to the fusion times detected
  in the \code{fusedanova} function.}
}
\description{
  Produce a prediction for a vector of \code{lambda}
  parameter and an array of \code{class}.
}
\examples{
data(aves)
fa <- fusedanova(x=aves$weight, class=aves$family, weight="laplace", gamma=5)
predict(fa, labels=aves$order)
}
\seealso{
  \code{\linkS4class{fusedanova}}.
}

