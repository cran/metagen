\name{formulaL}
\alias{formulaL}
\title{Regression coefficients: formulaL}
\usage{
  formulaL(y, d, h, g, x)
}
\arguments{
  \item{y}{k-vector of responses.}

  \item{d}{k-vector of heteroscedasticity.}

  \item{h}{scalar of heterogeneity.}

  \item{g}{p-vector of some p-variate Gaussian draw.}

  \item{x}{design k-p-matrix.}
}
\value{
  A p-vector.
}
\description{
  Calculate pivotal quantities for the regression
  coefficients using the method: formulaL form the
  dissertation.
}
\details{
  Algorithm for calculating a single generalised pivotal
  quantity for the regression coefficients for given
  generalised pivotal quantities for the heterogeneity
  using the univariate version of the pivotal formula.
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)

# When for example using the Mandel-Paule estimate:
bcg_h <- pfunc(y=bcg_y, d=bcg_d, x=bcg_x)(dim(bcg_x)[1] -
  dim(bcg_x)[2])

set.seed(51351) # for reproducibility
random_g <- rnorm(dim(bcg_x)[2])
formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg_x)

# The function can also be used when planing to perform
# a meta regression with no intercept, and only a singel
# covariate (i.e. dim(x) = 1).  In this case,
# the design matrix can simply be provided by a vector.
set.seed(51351) # for reproducibility
random_g <- rnorm(1)
formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=bcg$x)

# When performing a meta analysis, provide the function
# with a vector of 1s.
formulaL(y=bcg_y, d=bcg_d, h=bcg_h, g=random_g, x=rep(1,
  length(bcg_y)))
}

