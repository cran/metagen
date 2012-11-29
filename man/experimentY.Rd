\name{experimentY}
\alias{experimentY}
\title{Running a computer experiment}
\usage{
  experimentY(n, h, d, x, b, sgnf, piv_draws)
}
\arguments{
  \item{n}{number of draws.}

  \item{h}{heterogeneity.}

  \item{d}{heteroscedasticity.}

  \item{x}{design matrix.}

  \item{b}{regression coefficients.}

  \item{sgnf}{significance levels.}

  \item{piv_draws}{privotal draws.}
}
\value{
  Data frame of accumulated performance results.
}
\description{
  Runs a computer experiment that evaluates the performance
  of different inference methods for the random effects
  meta regression model with respect to heterogeneity and
  regression coefficients.
}
\examples{
h_test <- 0.03
x_test <- cbind(1,1:7)
b_test <- c(.5, .25)
sgnf_test <- c(0.025, 0.01)

set.seed(5133568) # for reproducibility
d_test <- rchisq(7, df=0.02)

# In an actual computer experiment, use 'piv_draws=1000' instead!!
experimentY(n=5, h=h_test, d=d_test, x=x_test, b=b_test,
  sgnf=sgnf_test, piv_draws=50)
}

