\name{rD}
\alias{rD}
\title{Data generation: Gaussian-Gaussian model}
\usage{
  rD(n, d, s)
}
\arguments{
  \item{n}{number of draws.}

  \item{d}{heteroscedasticity.}

  \item{s}{study sizes.}
}
\value{
  A (k,n)-matrix.  Each column is an independent draw.
}
\description{
  Random draws of heteroscedasticity responses of studies,
  where each study in a random effects meta regression
  model follows a Gaussian response. Thus D = (d * X) /
  (s-1) where X is chi-squared distributed.
}
\examples{
d_test = rchisq(13, df=0.02)
s_test = rep(100, 13)
rD(n=10, d=d_test, s=s_test)
}

