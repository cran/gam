\name{kyphosis}
\alias{kyphosis}
\docType{data}
\title{A classic example dataset for GAMs}
\description{
Data on the results of a spinal operation "laminectomy" on children, to
correct for a condition called "kyphosis"; see Hastie and
Tibshirani (1990) for details}
\usage{
data(kyphosis)
}
\format{A data frame with 81 observations on the following 4 variables.
  \describe{
    \item{Kyphosis}{a response factor with levels \code{absent}
      \code{present}.}
    \item{Age}{of child in months, a numeric vector}
    \item{Number}{of vertebra involved in the operation,a numeric vector}
    \item{Start}{level of the operation, a numeric vector}
  }
}
\source{
  Hastie, T. and Tibshirani, R. (1990)
  \emph{Generalized Additive Models.}
  London: Chapman and Hall.
}
\keyword{datasets}
