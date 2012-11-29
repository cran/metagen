\name{hConfidence}
\alias{hConfidence}
\title{Inference: Based on methods of moments and
maximum likelihood.}
\usage{
  hConfidence(y, d, x, sgnf)
}
\arguments{
  \item{y}{k-vector of study responses}

  \item{d}{k-vector of heteroscedasticity}

  \item{x}{design k-p-matrix}

  \item{sgnf}{significance levels}
}
\value{
  A data frame containing the bounds of the interval
  estimate.
}
\description{
  Calculates the so called Q-profiling confidence interval
  for the heterogeneity for data following a random effects
  meta regression model.
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_s <- bcg$size
bcg_x <- cbind(1,bcg$x)
sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)

set.seed(865287113) # for reproducibility

hConfidence(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025)
hConfidence(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev)
}

