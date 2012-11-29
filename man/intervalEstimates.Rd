\name{intervalEstimates}
\alias{intervalEstimates}
\title{Interval estimates: For the regression coefficients}
\usage{
  intervalEstimates(y, d, h_dat, x, sgnf)
}
\arguments{
  \item{y}{study responses.}

  \item{d}{heteroscedasticity.}

  \item{h_dat}{data frame of tau estimates.}

  \item{x}{design matrix.}

  \item{sgnf}{significance levels.}
}
\description{
  Interval estimates: For the regression coefficients
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)
bcg_h <- hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)
sgnf_lev <- c(0.01, 0.025, 0.05, 0.01)

intervalEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x, sgnf=0.025)
intervalEstimates(y=bcg_y, d=bcg_d, h_dat=bcg_h, x=bcg_x,
  sgnf=sgnf_lev)
}

