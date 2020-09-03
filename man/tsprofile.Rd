\name{tsprofile}
\alias{tsprofile}
\alias{tsprofile.tsvets.estimate}
\title{
Model Simulation Based Profiling
}
\description{
Profiling of model dynamics using simulation/estimation/prediction.
}
\usage{
\method{tsprofile}{tsvets.estimate}(object, h = 1, nsim = 100,
seed = NULL, cores = 1, trace = 0, solver = "nlminb", ...)
}
\arguments{
\item{object}{
An object of class tsvets.estimate
}
\item{h}{
The forecast horizon on which to evaluate performance metrics.
}
\item{nsim}{
The number of paths to generate.
}
\item{seed}{
An object specifying if and how the random number generator should be initialized.
See the simulate documentation for more details.
}
\item{cores}{
The number of cores to use for parallel estimation of the paths.
}
\item{trace}{
Whether to show a progress bar for the estimation of the paths.
}
\item{solver}{
Choice of solver to use for the estimation of the paths.
}
\item{\dots}{
Not currently used
}
}
\details{
The function profiles an estimated model by simulating and then estimating multiple paths from the
assumed DGP while leaving h values out for prediction evaluation. Each simulated path is equal to
the size of the original dataset plus h additional values, and initialized with the initial state
vector from the model. A data.table matrix is returned with the distribution of the coefficients 
from each path estimation as well as a stats table with the MAPE, BIAS and MSLRE by horizon, 
simulation and series.
}
\value{
See details above.
}
\author{
Alexios Galanos
}
