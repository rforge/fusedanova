\docType{methods}
\name{dataconvert}
\alias{dataconvert}
\alias{dataconvert,fusedanova-method}
\title{Conversion method for a fusedanova object}
\usage{
\\S4method{dataconvert}{fusedanova}(object, predicted=FALSE,
                  formating = c("df","list")[1], labels=FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{fusedanova}.}

  \item{predicted}{logical; if \code{TRUE}, the return
  value expands the \code{result} slot of the original
  \code{fusedanova} object. Otherwise, the
  \code{prediction} slot is explanded. Default is
  \code{FALSE}.}

  \item{formating}{a string, either \code{"list"} or
  \code{"df"}: does the return value should be a list of
  dataframes for each variable or only one dataframe with
  an additional column indexing the variables. By default,
  \code{"df"}.}

  \item{labels}{logical. Should the labels be sent back.
  Default is \code{FALSE}.}

  \item{...}{used for S4 compatibility.}
}
\description{
Convert the compressed fused-ANOVA path to a more handy but
more memory demanding data.frame format. This is typically
used for plotting purposes.
}
\examples{
data(aves)
fa <- fusedanova(x=aves$weight, class=aves$family)
dataconvert(fa)
}
\seealso{
\code{\linkS4class{fusedanova}}.
}

