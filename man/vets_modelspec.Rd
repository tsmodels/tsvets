\name{vets_modelspec}
\alias{vets_modelspec}
\title{
Model Specification
}
\description{
Specifies a vector ETS model prior to estimation.
}
\usage{
vets_modelspec(y, level = c("constant", "diagonal", "common", "full", 
"grouped"), slope = c("none", "constant", "common", "diagonal", 
"full", "grouped"), damped = c("none", "common", "diagonal", 
"full", "grouped"), seasonal = c("none", "common", "diagonal", 
"full", "grouped"), group = NULL, xreg = NULL, xreg_include = NULL, 
frequency = 1, lambda = NULL, lambda_lower = 0, lambda_upper = 1.5, 
dependence = c("diagonal", "full", "equicorrelation", "shrinkage"), 
cores = 1) 
}
\arguments{
\item{y}{
A multivariate xts matrix.
}
\item{level}{
The dynamics for the level component.
}
\item{slope}{
The dynamics for the slope component.
}
\item{damped}{
The dynamics for the dampening component.
}
\item{seasonal}{
THe dynamics for the seasonal component.
}
\item{group}{
A vector of indices denoting which group the series
belongs to (when using the grouped dynamics).
}
\item{xreg}{
An xts matrix of external regressors.
}
\item{xreg_include}{
A matrix of dimension ncol(y) by ncol(xreg) populated with
either 0, 1 or 2 (0 = no beta, 1 = individual beta and 
2 = grouped beta). 
}
\item{lambda}{
The Box Cox power transformation vector (see \code{box_cox}) in the
\bold{tsaux} package. If a single NA, then it will calculate optimal lambda 
based on the multivariate Box Cox approach of Velila (1993),
else if a vector of NA values it will calculate the individual Box Cox
optimal parameters. Can also be either a single value (common lambda) or
vector of values (individual lambda).
}
\item{lambda_lower}{
Lower bound on Box Cox lambda transformation.
}
\item{lambda_upper}{
Upper bound on Box Cox lambda transformation.
}
\item{frequency}{
The seasonal frequency.
}
\item{dependence}{
The dependence structure to impose.
}
\item{cores}{
The number of cores to use. This is used to initialize the states
of each series using the \code{tsets} package.
}
}
\details{
The specification allows to specify a vector additive damped ETS model
with options for the dynamics of the states and dependence.
}
\value{
An object of class \dQuote{tsvets.spec} with the following slots:\cr
\item{target}{A list with original data series, the data series index and the sampling
frequency}
\item{transform}{A list with details on the Box Cox transformation}
\item{model}{A list with details the type of model dynamics}
\item{dependence}{A list with details about the dependence structure}
\item{xreg}{A list with details on the external regressors}
\item{vets_env}{An environment with precalculated state matrices and other
parameters which will be passed to the estimation routine}
}
\references{
Athanasopoulos, G and de Silva, A. (2012),\emph{Multivariate Exponential Smoothing
for Forecasting Tourist Arrivals}, Journal of Travel Research 51(5) 640–-652.\cr
de Silva, A., R. Hyndman, and R. D. Snyder. (2010).\emph{The Vector
Innovations Structural Time Series Framework: A Simple Approach to Multivariate Forecasting},
Statistical Modelling (10) 353--74.
}
\author{
Alexios Galanos
}