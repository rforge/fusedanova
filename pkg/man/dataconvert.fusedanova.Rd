\docType{methods}
\name{dataconvert,fusedanova-method}
\alias{dataconvert.fusedanova}
\alias{dataconvert,fusedanova-method}
\title{Conversion method for a fusedanova object}
\arguments{
  \item{object}{an object of class \code{fusedanova}.}

  \item{predicted}{logical; if \code{TRUE}, the return is
  calculated on the \code{result} slot, else on the
  \code{prediction} slot. By default, \code{FALSE}.}

  \item{listformat}{logical; does the return should be a
  list of dataframes for each variable or only one
  dataframe with one more column indexing the variables. By
  default, \code{FALSE}.}

  \item{...}{used for S4 compatibility.}
}
\description{
  Produce for each object (dataframe, vector of oeder,
  vector of class) in the list of \code{result} or
  \code{prediction} (if \code{predicted} == \code{TRUE}) a
  new dataframe containing for each fusion the
  \code{lambda}, \code{beta}, \code{slope}, \code{class}.
}
\examples{
data(aves)
fa <- fusedanova(x=aves$weight, class=aves$family)
dataconvert(fa)
}
\seealso{
  \code{\linkS4class{fusedanova}}.
}

