% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.gnsstsmodel}
\alias{print.gnsstsmodel}
\title{Print method for a \code{gnsstsmodel} object.}
\usage{
\method{print}{gnsstsmodel}(x, ...)
}
\arguments{
\item{x}{A \code{gnsstsmodel} object.}

\item{...}{Additional graphical parameters.}
}
\value{
No return value. Print a \code{gnsstsmodel} object.
}
\description{
Print method for a \code{gnsstsmodel} object.
}
\examples{
\dontrun{
data(cola)
fit_gmwmx = estimate_gmwmx(x = cola,
                           theta_0 = c(0.1,0.1,0.1,0.1), 
                           n_seasonal = 1, 
                           model_string = "wn+matern")
                           
print(fit_gmwmx)
}

}
