\name{plot}
\alias{plot}
\alias{plot.tsvets.estimate}
\alias{plot.tsvets.predict}
\title{
Estimated Model Plots
}
\description{
Plots for objects generated from the tsvets functions.
}
\usage{
\method{plot}{tsvets.estimate}(x, y = NULL, 
type = c("fitted", "states", "residuals"), 
series = 1:min(10, ncol(fitted(x))), ...)
\method{plot}{tsvets.predict}(x, y = NULL, 
series = 1, n_original = NULL, ...) 
}
\arguments{
\item{x}{
An object of class \dQuote{tsvets.estimate}.
}
\item{y}{
Not used.
}
\item{type}{
The type of plot with a choice of \dQuote{fitted},
\dQuote{states} and \dQuote{residuals}.
}
\item{series}{
The selection of series to plot, with a maximum of 10
at a time for object of class \dQuote{tsvets.estimate},
and the single series for the object of class
\dQuote{tsvets.predict}
}
\item{n_original}{
Number of actual data points to include in the predictio 
plot.
}
\item{\dots}{
Additional arguments passed to the underlying plot function.
}
}
\author{
Alexios Galanos
}
