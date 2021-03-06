\name{tsbenchmark}
\alias{tsbenchmark}
\alias{tsbenchmark.tsvets.spec}
\title{
Model Benchmarking
}
\description{
Benchmarks and saves the details of a model for benchmarking of speed
and unit testing.
}
\usage{
\method{tsbenchmark}{tsvets.spec}(object, solver = "nlminb", 
control = list(trace = 0), ...)
}
\arguments{
\item{object}{
An object of class tsvets.spec
}
\item{solver}{
The solver to use for estimation.
}
\item{control}{
Any additional controls passed to the solver.
}
\item{\dots}{
Not currently used
}
}
\details{
This is meant mostly for internal development use.
}
\value{
A data.table
}
\author{
Alexios Galanos
}
