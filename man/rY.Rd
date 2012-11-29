\name{rY}
\alias{rY}
\title{Data generation: Gaussian-Gaussian model}
\usage{
  rY(n, h, d, x, b)
}
\arguments{
  \item{n}{number of draws.}

  \item{h}{heterogeneity.}

  \item{d}{heteroscedasticity.}

  \item{x}{design matrix.}

  \item{b}{regression coefficients.}
}
\value{
  A (k,n)-matrix.  Each column is an independent draw.
}
\description{
  Random draws of response vectors y following the
  distribution of a random effects meta regression model.
  Each column is an independent draw.
}
\examples{
x_test = cbind(1,1:13)
h_test = .03
d_test = rchisq(13, df=0.02)
b_test = c(0.02, 0.03)
rY(n=10, h=h_test, d=d_test, x=x_test, b=b_test)
}

