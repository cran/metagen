\name{regressionEstimates}
\alias{regressionEstimates}
\title{Point estimates: For the regression coefficients}
\usage{
  regressionEstimates(y, d, h_dat, x)
}
\arguments{
  \item{y}{study responses, k-vector of responses.}

  \item{d}{heteroscedasticity, k-vector of
  heteroscedasticities.}

  \item{h_dat}{Here, 'h_dat' should be a data frame with
  variables 'type' and 'h'.  Thus, one may use h_dat =
  hEstimates(y, d, x).}

  \item{x}{design matrix, k-p-matrix.}
}
\value{
  A list of estimates for the regression coefficients.

  Here, 'h_dat' should be a data frame with variables
  'type' and 'h', thus, we may use h_dat = hEstimates(y, d,
  x)
}
\description{
  Calculates point estimates for the regression coefficient
  for given point estimates of the variance components 'd'
  and a data frame of different estimates of the
  heterogeneity 'h'.
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)
bcg_h <- hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)
regressionEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x)
}

