\name{rB}
\alias{rB}
\title{Data generation: Log-risk-ration of a binomial-Gaussian model}
\usage{
  rB(n, h, s, a, r, x, b)
}
\arguments{
  \item{n}{number of draws.}

  \item{h}{heterogeneity.}

  \item{s}{study sizes.}

  \item{a}{balance of group assignments.}

  \item{r}{fixed risk in the treatment group.}

  \item{x}{design matrix.}

  \item{b}{regression coefficients.}
}
\value{
  A (2k,n) matrix.  Each column is an independent draw.
}
\description{
  Random draws of log risk ratios from a hierarchical
  binomial Gaussian model.
}
\details{
  It is always assumed that at least one response in a
  study has happend, i.e., a response of 0 in a treatment
  or control group is rounded up to 1.  Note that this may
  lead to an overestimation of small risks.  If possible,
  make sure your sample sizes are large enough to
  compensate for this effect.
}
\examples{
h_test <- .03
x_test <- cbind(1,1:13)
b_test <- c(0.02, 0.03)
s_test <- rep(2000, 13)
a_test <- rep(.3, 13)
rB(n=10, h=h_test, s=s_test, a=a_test, r=.3, x=x_test, b=b_test)
}

