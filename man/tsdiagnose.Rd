\name{tsdiagnose}
\alias{tsdiagnose.tsvets.estimate}
\alias{tsdiagnose}
\title{
Model Diagnostics
}
\description{
Creates a short summary of model based diagnosticts.
}
\usage{
\method{tsdiagnose}{tsvets.estimate}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{\dots}{
Not currently used.
}
}
\value{
A list of tables (printed out and returned insivibly) with
diagnostics including the eigenvalues of the D matrix as
well as Normality and multivariate outlier tests.
}
\author{
Alexios Galanos.
}
