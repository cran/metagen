\name{experimentD}
\alias{experimentD}
\title{Running a computer experiment}
\usage{
  experimentD(n, h, d, s, x, b, sgnf, piv_draws)
}
\arguments{
  \item{n}{number of draws.}

  \item{h}{heterogeneity.}

  \item{d}{heteroscedasticity.}

  \item{s}{vector study sizes.}

  \item{x}{design matrix.}

  \item{b}{regression coefficients.}

  \item{sgnf}{significance levels.}

  \item{piv_draws}{privotal draws.}
}
\value{
  Data frame of accumulated analysis of random data drawn
  independend from the respeced random effects meta
  regression model.
}
\description{
  Runs a computer experiment that evaluates the performance
  of different inference methods for the random effects
  meta regression model with respect to heterogeneity and
  regression coefficients.
}
\details{
  This also includes methods adjusting for uncertainty in
  the heteroscedasticity vector.  In particular, the study
  sizes need to be known, here.
}
\examples{
h_test <- 0.03
x_test <- cbind(1,1:7)
b_test <- c(.5, .25)
sgnf_test <- c(0.025, 0.01)

set.seed(5133568) # for reproducibility
d_test <- rchisq(7, df=0.02)
s_test <- runif(7, min=200, max=2000)

# In an actual computer experiment, use 'piv_draws=1000' instead!!
experimentD(n=5, h=h_test, d=d_test, s=s_test, x=x_test, b=b_test,
  sgnf=sgnf_test, piv_draws=50)
}

