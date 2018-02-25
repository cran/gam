\name{gam.exact}
\alias{gam.exact}
\title{A method for gam producing asymptotically exact standard errors
  for linear estimates}
\description{This function is a "wrapper" for a Gam object, and produces exact standard errors
       for each linear term in the gam call (except for the intercept).}
\usage{
gam.exact(Gam.obj)
}
\arguments{
  \item{Gam.obj}{a Gam object}
}
\details{
 Only standard errors for the linear terms are produced.
There is a print method for the Gamex class.
}
\value{
 A list (of class Gamex) containing a table of coefficients and a variance
 covariance matrix for the linear terms in the formula of the gam call.
 }
\references{[1] Issues in Semiparametric Regression: A Case Study of Time Series Models
           in Air Pollution and Mortality,  Dominici F., McDermott A., Hastie T.J.,
           \emph{JASA}, December 2004, 99(468), 938-948. See
	   \url{http://web.stanford.edu/~hastie/Papers/dominiciR2.pdf}
	 }
\author{Aidan McDermott, Department of Biostatistics, Johns
	   Hopkins University. Modified by Trevor Hastie for R}

\examples{
set.seed(31)
n     <- 200
x     <- rnorm(n)
y     <- rnorm(n)
a     <- rep(1:10,length=n)
b     <- rnorm(n)
z     <- 1.4 + 2.1*a + 1.2*b + 0.2*sin(x/(3*max(x))) + 0.3*cos(y/(5*max(y))) + 0.5 * rnorm(n)
dat   <- data.frame(x,y,a,b,z,testit=b*2)
### Model 1: Basic
Gam.o <- gam(z ~ a + b + s(x,3) + s(y,5), data=dat)
coefficients(summary.glm(Gam.o))
gam.exact(Gam.o)
### Model 2: Poisson
Gam.o <- gam(round(abs(z)) ~ a + b + s(x,3) + s(y,5), data=dat,family=poisson)
coefficients(summary.glm(Gam.o))
gam.exact(Gam.o)
}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
