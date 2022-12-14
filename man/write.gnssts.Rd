% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnssts.R
\name{write.gnssts}
\alias{write.gnssts}
\title{Write a gnssts object}
\usage{
write.gnssts(x, filename, format = "mom")
}
\arguments{
\item{x}{A \code{R} object to save as a \code{gnssts} object.}

\item{filename}{A \code{string} specifying the name of the file to write.}

\item{format}{A \code{string} specifying the format of the file to write.}
}
\value{
No return value. Write a \code{gnssts} object in a .mom file by default.
}
\description{
Write a gnssts object
}
\examples{
phase <- 0.45
amplitude <- 2.5
sigma2_wn <- 15
bias <- 0
trend <- 5 / 365.25
cosU <- amplitude * cos(phase)
sinU <- amplitude * sin(phase)
year <- 5
n <- year * 365
jump_vec <- c(200, 300, 500)
jump_height <- c(10, 15, 20)
nbr_sin <- 1
A <- create_A_matrix(1:n, jump_vec, n_seasonal = nbr_sin)
x_0 <- c(bias, trend, cosU, sinU, jump_height)
eps <- rnorm(n = n, sd = sqrt(sigma2_wn))
yy <- A \%*\% x_0 + eps
gnssts_obj <- create.gnssts(t = 1:length(yy), y = yy, jumps = jump_vec)
str(gnssts_obj)
\dontrun{
write.gnssts(x = gnssts_obj, filename = "test.mom")
}
}
