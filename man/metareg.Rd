\name{metareg}
\alias{metareg}
\title{Inference: Based on methods of moments and
maximum likelihood.}
\usage{
  metareg(y, d, x, sgnf)
}
\arguments{
  \item{y}{k-vector of study responses}

  \item{d}{k-vector of heteroscedasticity}

  \item{x}{design k-p-matrix}

  \item{sgnf}{significance levels}
}
\value{
  The same return type as the skeleton 'metagenEmpty()'.
}
\description{
  Calculates common statistics for point and confidence
  interval estimates for the heterogeneity and the
  regression coefficients of the random effects meta
  regression model based on the given data.
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_s <- bcg$size
bcg_x <- cbind(1,bcg$x)
sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)

set.seed(865287113) # for reproducibility

c1 <- metareg(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=0.025)
c2 <- metareg(y=bcg_y, d=bcg_d, x=bcg_x, sgnf=sgnf_lev)

if (!all(names(c1) == names(metagenEmpty()))) stop("Name clash")
if (!all(names(c2) == names(metagenEmpty()))) stop("Name clash")
}

