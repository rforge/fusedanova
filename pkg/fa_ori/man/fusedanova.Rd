\name{fusedanova}
\alias{fusedanova}
\title{Fit a Fused ANOVA model}
\usage{
fusedanova(x, class, ...)
}
\arguments{
  \item{x}{matrix (or column vector) which rows represent
  individuals and columns independant variables.}

  \item{class}{vector or factor giving the class of each
  individual. If missing, \code{1:nrow(x)} is used
  (clustering mode with one individual per class).}

  \item{...}{list of additional parameters to overwrite the
  defaults of the fitting procedure. Include : \itemize{%

  \item{\code{weights}: } {character; which type of weights
  is supposed to be used.  The supported weights are :
  \code{"default"}, \code{"laplace"}, \code{"gaussian"},
  \code{"adaptive"}, \code{"naivettest"}, \code{"ttest"},
  \code{"welch"} and \code{"personal"}. Se details below.
  By default, its value is \code{"default"}.}

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

The optimization problem solved by fused-ANOVA is
\if{latex}{\deqn{% \hat{\beta}_{\lambda} = \arg
\min_{\beta} \left\{\sum_{k=1}^K \sum_{i=1}^{n_k}
\left(Y_{ik}-\beta_k \right)^2 + \lambda \sum_{k,\ell}
w_{kl} \left|\beta_k - \beta_\ell \right|\right\}}}
\if{html}{\out{ <center> &beta;<sup>hat</sup>
<sub>&lambda;<sub>1</sub></sub> = argmin<sub>&beta;</sub>
sum<sub>k</sub> sum_i (Y<sub>ik</sub> -
&beta<sub>k</sub>)<sup>2</sup> + &lambda; sum<sub>k,l</sub>
w<sub>k,l</sub> &#124; &beta;<sub>k</sub> -
&beta;<sub>l</sub> &#124;, </center> }}
\if{text}{\deqn{beta.hat(lambda) = argmin_beta sum_k sum_i
(Y_ik - beta_k)^2 + lambda sum_k sum_l w_kl | beta_k -
beta_l|,}}

where \eqn{Y_{ik}}{Y_ik} is the intensity of a continuous
random variable for sample \eqn{i}{i} in condition
\eqn{k}{k} and \eqn{\beta_k}{beta_k} is the mean parameter
of condition \eqn{k}{k}. We denote by \eqn{K}{K} the total
number of conditions and \eqn{n_k}{n_k} the number of
sample in each condition.

More details related to the weights are coming...
}
\description{
Adjust a penalized ANOVA model with Fused-LASSO (or Total
Variation) penality, ie. a sum of weighted
\eqn{\ell_1}{l1}-norm on the difference of each
coefficient. See details below.
}
\examples{
\dontrun{
data(aves)
fa.laplace <- fusedanova(x=aves$weight, class=aves$family, weights="laplace", gamma=5)
plot(fa.laplace, labels=aves$order)

fa.ttest <- fusedanova(x=aves$weight, class=aves$family, weights="naivettest")
plot(fa.ttest, labels=aves$order)

fa.ada <- fusedanova(x=aves$weight, class=aves$family, weights="adaptive", gamma=2)
plot(fa.ada, labels=aves$order)
}
}
\seealso{
See also \code{\linkS4class{fusedanova}},
\code{\link{plot,fusedanova-method}} and
\code{\link{cv.fa}}.
}
\keyword{models,}
\keyword{regression}

