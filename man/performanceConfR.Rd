\name{performanceConfR}
\alias{performanceConfR}
\title{Running a computer experiment: Adding performance measures}
\usage{
  performanceConfR(accum_int, true)
}
\arguments{
  \item{accum_int}{accumulated interval estimates.  At
  least the following columns need to be present: lower and
  upper and parameter.}

  \item{true}{true parameter.}
}
\description{
  Adding performance measurements to accumulated results of
  a computer experiment running multiple analysis of
  different simulated data following a random effects meta
  regression model.
}
\details{
  Adds performance measurements to interval estimates of
  the regression coefficients.
}
\examples{
# For an example, see the 'performance' function.
}

