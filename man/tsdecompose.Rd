\name{tsdecompose}
\alias{tsdecompose}
\alias{tsdecompose.tsvets.estimate}
\title{
Model Decomposition
}
\description{
Decomposes the estimated model or prediction into its component parts.
}
\usage{
\method{tsdecompose}{tsvets.estimate}(object, simplify = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{simplify}{
Simplification of the returned states aggregates the level and slope 
(if present) into a Trend, Seasonal, X and Irregular components.
}
\item{\dots}{
Not currently used
}
}
\value{
Returns a list of xts matrices with the state components (including error).
}
\author{
Alexios Galanos
}
