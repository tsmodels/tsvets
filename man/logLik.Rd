% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{logLik.tsvets.estimate}
\alias{logLik.tsvets.estimate}
\alias{logLik}
\title{Model Log-Likelihood}
\usage{
\method{logLik}{tsvets.estimate}(object, ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsvets.estimate}.}

\item{...}{not currently used.}
}
\value{
Returns an object of class logLik. This is a number with at least one
attribute, "df" (degrees of freedom), giving the number of (estimated)
parameters in the model.
}
\description{
Extract the log-likelihood from an estimated model.
}
