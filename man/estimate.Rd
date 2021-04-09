\name{estimate}
\alias{estimate.tsvets.spec}
\title{Model Estimation}
\description{
Estimates a model given a specification object using maximum likelihood.
}
\usage{
\method{estimate}{tsvets.spec}(object, solver = "nlminb", 
control = list(trace = 1, iter.max = 200, eval.max = 1000), ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.spec}.
}
\item{solver}{
Currently \dQuote{L-BFGS-B} from \dQuote{optim}, \dQuote{nlminb}, \dQuote{solnp}, \dQuote{gosolnp} and \dQuote{nloptr} 
bound constrained solvers are supported
}
\item{control}{
Control parameters for the solvers.
}
\item{\dots}{
Currently allowed for passing additional arguments to the \dQuote{gosolnp} solver.
}
}
\details{
Minimization of the negative of the log likelihood function of the vector Additive
ETS model with a soft barrier based on the stability constraint of the models.
}
\value{
An object of class \dQuote{tsvets.estimate}.
}
\references{
Athanasopoulos, G and de Silva, A. (2012),\emph{Multivariate Exponential Smoothing
for Forecasting Tourist Arrivals}, Journal of Travel Research 51(5) 640–-652.\cr
de Silva, A., R. Hyndman, and R. D. Snyder. (2010).\emph{The Vector
Innovations Structural Time Series Framework: A Simple Approach to Multivariate Forecasting},
Statistical Modelling (10) 353--74.
}
\author{
Alexios Galanos
}
