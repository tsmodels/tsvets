\name{fitted}
\alias{fitted.tsvets.estimate}
\title{
Model Fitted Values
}
\description{
Extract the fitted values from an estimated model.
}
\usage{
\method{fitted}{tsvets.estimate}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{\dots}{
Not currently used.
}
}
\details{
An xts matrix, back adjusted if a transformation was used.
}
\value{
An xts matrix.
}
\author{
Alexios Galanos
}
