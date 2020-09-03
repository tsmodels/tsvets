\name{residuals}
\alias{residuals}
\alias{residuals.tsvets.estimate}
\title{
Model Residuals
}
\description{
Extract the residual values from an estimated model.
}
\usage{
\method{residuals}{tsvets.estimate}(object, raw = FALSE, ...)
}
\arguments{
\item{object}{
An object of class tsvets.estimate
}
\item{raw}{
If TRUE returns the untransformed residuals.
}
\item{\dots}{
Not currently used
}
}
\value{
An xts vector of the model residuals.
}
\author{
Alexios Galanos
}
