\name{tsbacktest}
\alias{tsbacktest}
\alias{tsbacktest.tsvets.spec}
\title{
Walk Forward Model Backtest
}
\description{
Generates an expanding window walk forward backtest.
}
\usage{
\method{tsbacktest}{tsvets.spec}(object, start = floor(NROW(object$target$y_orig)/2), 
end = NROW(object$target$y_orig), h = 1, alpha = NULL, 
aggregate = FALSE, weights = NULL, 
cores = 1, save_output = FALSE, save_dir = "~/tmp/", 
solver = "nlminb", trace = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsvets.spec}.
}
\item{start}{
The numeric data index from which to start the backtest.
}
\item{end}{
The numeric data index on which to end the backtest. The backtest
will end 1 period before that date in order to have at least 1 out
of sample value to compare against.
}
\item{h}{
The forecast horizon. As the expanding window approaches the \dQuote{end},
the horizon will automatically shrink to the number of available out of
sample periods.
}
\item{alpha}{
An optional numeric vector of coverage rates for which to calculate
the quantiles.
}
\item{aggregate}{
Whether to create an aggregate series from the forecasts (see 
\code{\link{tsaggregate}} for details)
}
\item{weights}{
A numeric weights vector for the aggregation.
}
\item{cores}{
The number of cores to use for parallel estimation. The package
makes use of the doSNOW package.
}
\item{save_output}{
Whether to save each estimated and predicted object.
}
\item{save_dir}{
The directory to save the estimated and predicted objects.
These are saved as rds files using \dQuote{model_estimationDate}
and \dQuote{predict_estimationDate} as names with estimationDate
the date on which the model was estimated.
}
\item{solver}{
The solver to use.
}
\item{trace}{
Whether to show the progress bar.
}
\item{\dots}{
Not currently used
}
}
\value{
A list with the following data.tables:
\itemize{
\item prediction : the backtest table with forecasts and actuals by series
\item metrics: a summary performance table showing metrics by 
forecast horizon and series
}
}
\details{
Note that unless NULL is specified for lambda in the object, it
will be estimated for each estimation period. More flexibility
may be provided in future releases.
}
\author{
Alexios Galanos
}
