\name{tsmetrics}
\alias{tsmetrics.tsvets.estimate}
\alias{tsmetrics.tsvets.predict}
\title{
Model Evaluation Metrics
}
\description{
Generates a set of model evaluation metrics.
}
\usage{
\method{tsmetrics}{tsvets.estimate}(object, weights = NULL, ...)
\method{tsmetrics}{tsvets.predict}(object, actual, weights = NULL, alpha = 0.1, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{weights}{
A vector of series importance weights (which sum to 1), used
in generating weighted performance metrics.
}
\item{actual}{
The actual data matched to the dates of the forecasts.
}
\item{alpha}{
The coverage level for distributional forecast metrics.
}
\item{\dots}{
Not currently used.
}
}
\value{
A one row data.frame containing summary performance metrics.
}
\author{
Alexios Galanos
}
