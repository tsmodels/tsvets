\name{tsfilter.tsvets.estimate}
\alias{tsfilter}
\alias{tsfilter.tsvets.estimate}
\title{
Online Model Filtering
}
\description{
Online filter which updates the states and fitted values using new data.
}
\usage{
\method{tsfilter}{tsvets.estimate}(object, y, newxreg = NULL, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{y}{
An xts matrix of new information related to y.
The function checks whether y contains indices (dates)
which are not in the passed object and only filters
new information.
}
\item{newxreg}{
An xts matrix of new information related to external
regressors (if those were used in the original model
estimated).
}
\item{\dots}{
Not currently used.
}
}
\details{
The new data is filtered (1 step ahead) using the last
state of the object. Once this is complete, the object
is updated with the new states and information so that
the process can be continued on the same object as new
information arrives.
}
\value{
An object of class \dQuote{tsvets.estimate}
}
\author{
Alexios Galanos
}
