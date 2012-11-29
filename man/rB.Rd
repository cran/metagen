\name{rB}
\alias{rB}
\title{Data generation: Log-risk-ration of a Binomial-Gaussian model}
\usage{
  rB(n, h, s, a, r, x, b)
}
\arguments{
  \item{n}{number of draws}

  \item{h}{heterogeneity}

  \item{s}{study sizes}

  \item{a}{balance of group assignments}

  \item{r}{fixed risk in the treatment group}

  \item{x}{design matrix}

  \item{b}{regression coefficients}
}
\value{
  A (2k,n) matrix.  Each column is an independent draw.
  x_test <- cbind(1,1:13) h_test <- .03 a_test <- runif(13,
  min-.3, max=.3) b_test <- c(0.02, 0.03) rB(n=10,
  h=h_test, a=a_test, r=.3, x=x_test, b=b_test)
}
\description{
  Random draws of log risk ratios from a hierarchical
  Binomial Gaussian model.
}
\details{
  Same comment applies as in the function 'rBinomGauss'.
}

