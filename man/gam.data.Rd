\name{gam.data}
\alias{gam.data}
\alias{gam.newdata}
\docType{data}
\title{Simulated dataset for gam}
\description{
A simple simulated dataset, used to test out the gam functions
}
\usage{
data(gam.data)
data(gam.newdata)
}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{x}{a numeric vector - predictor}
    \item{y}{a numeric vector - the response}
    \item{z}{a numeric vector - noise predictor}
    \item{f}{a numeric vector - true function}
    \item{probf}{a numeric vector - probability function}
    \item{ybin}{a numeric vector - binary response}
  }
}
\details{
This dataset is artificial, and is used to test out some of the features
of gam.
}
\examples{
data(gam.data)
gam(y ~ s(x) + z, data=gam.data)
}
\keyword{datasets}
