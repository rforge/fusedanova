\docType{class}
\name{fusedanova-class}
\alias{fusedanova-class}
\alias{predict,fusedanovaquadrupen-method}
\alias{print,fusedanova-method}
\alias{show,fusedanova-method}
\title{Class "fusedanova"}
\description{
  Class of object returned by the fusedanova function.
}
\section{Slots}{
  \describe{

  \item{\code{result}: }{a list of objects, one per
  variable. Objects are lists containing : \itemize{%
  \item{\code{table}: }{a dataframe containing the
  \code{lambda}; \code{Beta}, value of coefficient ;
  \code{slope}; \code{idown} the minimal index of class in
  the group and \code{iup}, the maximum one; for each
  fusion time \code{lambda}.}

  \item{\code{order}: }{a vector giving the order of means
  before any fuse. Needed with \code{classes} for data
  conversion to a more understandable format.}

  } }

  \item{\code{weights}: }{the weights used t eprform the
  fit.}

  ' \item{\code{classes}: }{the intitulate of the classes
  as entered in the \code{class} parameter in the
  \code{fusedanova} function.}

  \item{\code{prediction}: }{if a \code{lambdalist} was
  given to the \code{fusedanova} function, a list contaning
  for each variable an object with the same format as the
  ones in \code{result}. Else, an empty list.}

  \item{\code{lambdalist}: }{the \code{lambdalist} vector
  parameter given in the \code{fusedanova} function. If not
  given, empty vector.}

  \item{\code{algorithm}: }{a character string indicating
  wether splits where allowed or not.}

  }
}

\section{Methods}{
  Specific plotting, predict and conversion methods are
  available and documented (\code{\link{plot.fusedanova}},
  \code{\link{predict.fusedanova}},
  \code{\link{dataconvert.fusedanova}}).
}
\seealso{
  See also \code{\link{plot.fusedanova}},
  \code{\link{predict.fusedanova}}
  \code{\link{dataconvert.fusedanova}} and
  \code{\link{fusedanova}}.
}
\keyword{class}

