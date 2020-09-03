\name{tsaggregate}
\alias{tsaggregate}
\alias{tsaggregate.tsvets.estimate}
\alias{tsaggregate.tsvets.predict}
\alias{tsaggregate.tsvets.simulate}
\title{
Series Aggregation
}
\description{
Aggregates estimated, predicted or simulated series based on a weighting vector.
}
\usage{
\method{tsaggregate}{tsvets.estimate}(object, weights = NULL, return_model = FALSE, ...) 
\method{tsaggregate}{tsvets.predict}(object, weights = NULL, ...)
\method{tsaggregate}{tsvets.simulate}(object, weights = NULL, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}. \dQuote{tsvets.predict} or \dQuote{tsvets.simulate}.
}
\item{weights}{
A vector of weights of length equal to the number of series.
}
\item{return_model}{
If the estimated object is a homogeneous cofficients model with common lamda parameter (of either
0, 1 or NULL), then it will return a univariate tsvets object of class \dQuote{tsvets.estimate}.
}
\item{\dots}{
Not currently used
}
}
\value{
Depends on the input class.
}
\details{
For an estimated object which has common components (homogeneous coefficients) for all 
estimated states and common lambda parameter, then a reconstructed object of 
the same class representing the weighted representation of the model is returned.
In all other cases and input classes, the returned object will depend on whether
the lambda parameter was 0 or 1 for all underlying series. For the case of the
log transform (lambda = 0), then the states aggregate (given a weighting vector) 
whilst the actual, fitted, predicted or simulated values aggregate in logs and 
are then exponentiated (a multiplicative model). In all other cases only the actual, 
fitted, predicted or simulated values are returned representing the weighted 
aggregation of the underlying series given the weighting vector.
}
\author{
Alexios Galanos
}
\references{
Hyndman, Rob and Koehler, Anne B and Ord, J Keith and Snyder, Ralph D, 2008,
Forecasting with exponential smoothing: the state space approach, Section 17.1.2,
,Springer Science \& Business Media.
}