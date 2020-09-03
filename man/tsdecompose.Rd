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
\method{tsdecompose}{tsvets.estimate}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{\dots}{
Not currently used
}
}
\value{
Returns a list of xts matrices with the fitted and state components (including error).
}
\author{
Alexios Galanos
}
