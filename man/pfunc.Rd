\name{pfunc}
\alias{pfunc}
\title{The p_delta(eta) function.}
\usage{
  pfunc(y, d, x)
}
\arguments{
  \item{y}{study responses.}

  \item{d}{heteroscedasticity.}

  \item{x}{design matrix.}
}
\value{
  A vector valued function.
}
\description{
  Returns the p-function. (fastest version tested)
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)
pfunc(y=bcg_y, d=bcg_d, x=bcg_x)

# Calculating the Mandel-Paule estimate:
pfunc(y=bcg_y, d=bcg_d, x=bcg_x)(dim(bcg_x)[1] - dim(bcg_x)[2])
}

