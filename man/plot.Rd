% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.tsvets.estimate}
\alias{plot.tsvets.estimate}
\alias{plot}
\alias{plot.tsvets.predict}
\title{Object Plots}
\usage{
\method{plot}{tsvets.estimate}(
  x,
  y = NULL,
  type = c("fitted", "states", "residuals"),
  series = 1:min(10, ncol(fitted(x))),
  ...
)

\method{plot}{tsvets.predict}(x, y = NULL, series = 1, n_original = NULL, ...)
}
\arguments{
\item{x}{an object of class \dQuote{tsvets.estimate} or \dQuote{tsvets.predict}.}

\item{y}{not used.}

\item{type}{the type of plot with a choice of \dQuote{fitted}, \dQuote{states} 
and \dQuote{residuals}.}

\item{series}{the selection of series to plot, with a maximum of 10 at a time 
for object of class \dQuote{tsvets.estimate}, and the single series for the 
object of class \dQuote{tsvets.predict}}

\item{...}{additional arguments passed to the underlying plot function.}

\item{n_original}{number of actual data points to include in the prediction 
plot.}
}
\description{
Plots for objects generated from the tsvets functions.
}
