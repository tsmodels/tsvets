\name{summary}
\alias{summary}
\alias{summary.tsvets.estimate}
\title{
Model Estimation Summary
}
\description{
Summary method for class \dQuote{tsvets.estimate}
}
\usage{
\method{summary}{tsvets.estimate}(object, weights = NULL, ...)
}
\arguments{
\item{object}{
An object of class tsvets.estimate
}
\item{weights}{
A vector of weights (summing to 1) which is passed to the \code{tsmetrics}
method for generating weighted performance metrics.
}
\item{\dots}{
Not currently used
}
}
\value{
A printout of the parameter summary, model type and some model metrics.
}
\author{
Alexios Galanos
}
