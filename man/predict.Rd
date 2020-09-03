\name{predict}
\alias{predict.tsvets.estimate}
\title{
Model Prediction
}
\description{
Prediction function for class \dQuote{tsvets.estimate}.
}
\usage{
\method{predict}{tsvets.estimate}(object, h = 12, newxreg = NULL, nsim = 1000, forc_dates = NULL, 
init_states = NULL, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.estimate}.
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{nsim}{
The number of simulations to use for generating the simulated
predictive distribution.
}
\item{forc_dates}{
An optional vector of forecast dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular
sequence of dates after the last available date in the data.
}
\item{init_states}{
An optional vector of states to initialize the forecast. If NULL,
will use the last available state from the estimated model.
}
\item{\dots}{
No currently used.
}
}
\details{
Generates the predicted distribution of the vector Additive ETS model
simulating from the multivariate Normal distribution with covariance
based on the estimated model residuals.
}
\value{
An object of class \dQuote{tsvets.predict} with most relevant output a data.table which 
contains one row for each series and columns of tsmodel.predict objects, one 
for the predicted values, and for each of the states.
}
\author{
Alexios Galanos
}
