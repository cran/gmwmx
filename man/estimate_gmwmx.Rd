% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmwmx.R
\name{estimate_gmwmx}
\alias{estimate_gmwmx}
\title{Estimate a stochastic model in a two-steps procedure using the GMWMX estimator.}
\usage{
estimate_gmwmx(
  x,
  theta_0,
  n_seasonal = 1,
  model_string,
  method = "L-BFGS-B",
  maxit = 1e+06,
  ci = FALSE,
  k_iter = 1
)
}
\arguments{
\item{x}{A \code{gnssts} object}

\item{theta_0}{A \code{vector} specifying the initial values for the vector of parameter of the stochastic model considered.}

\item{n_seasonal}{An \code{integer} specifying the number of seasonal component in the time series.}

\item{model_string}{A \code{string} specifying the model to be estimated.}

\item{method}{A \code{string} specifying the numerical optimization method that should be supplied to \code{optim()}}

\item{maxit}{An \code{integer} specifying the maximum number of iterations for the numerical optimization procedure.}

\item{ci}{A \code{boolean} specifying if confidence intervals for the estimated parameters should be computed.}

\item{k_iter}{An \code{integer} specifying the number of time the two steps GMWMX procedure should be run.}
}
\value{
A \code{gnsstsmodel} object.
}
\description{
Estimate a stochastic model in a two-steps procedure using the GMWMX estimator.
}
\examples{
\dontrun{
data(cola)
fit_gmwmx = estimate_gmwmx(x = cola,
                           theta_0 = c(0.1,0.1,0.1,0.1), 
                           n_seasonal = 1, 
                           model_string = "wn+matern")
}

}
