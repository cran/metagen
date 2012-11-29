\name{rBinomGauss}
\alias{rBinomGauss}
\title{Data generation: Sampling data of clinical trials}
\usage{
  rBinomGauss(h, s, a, r, x, b)
}
\arguments{
  \item{h}{heterogeneity}

  \item{s}{study sizes}

  \item{a}{balance of group assignments}

  \item{r}{fixed risk in the control group}

  \item{x}{design matrix}

  \item{b}{regression coefficients}
}
\value{
  A list containing the risk and a data frame with the
  studies.
}
\description{
  A random draw of a hierarchical Binomial Gaussian model.
}
\details{
  We always assume at least one response in a study has
  happend, i.e. a response of 0 in a treatment or control
  group is rounded up to 1.  Note that this may lead to an
  overestimation of small risks.

  You may work around this by increasing study sizes.
}
\examples{
h_test <- .03
x_test <- cbind(1,1:13)
b_test <- c(0.02, 0.03)
s_test <- rep(2000, 13)
a_test <- rep(.3, 13)
dat <- rBinomGauss(  h=h_test, s=s_test, a=a_test, r=0.03
            , x=x_test, b=b_test)$study

if(!all(dim(dat) == c(dim(x_test)[1], 4))) stop("Wrong dimension")
}

