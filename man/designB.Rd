\name{designB}
\alias{designB}
\title{Design: Binomial responses}
\usage{
  designB(n, h_bounds, a_bounds, s_bounds, r, x)
}
\arguments{
  \item{n}{resolution of the heterogeneity.  n is the
  number of of different heterogeneity parameters in the
  design.}

  \item{h_bounds}{bounds of the heterogeneity.}

  \item{a_bounds}{bounds of the balancing factor of group
  assignments.}

  \item{s_bounds}{bounds of the study sizes.}

  \item{r}{fixed risk in the control.}

  \item{x}{design matrix.}
}
\value{
  Function returns a data frame.  Each line of this data
  frame can be an input to the function 'rB' which is used
  to sample data from such a design.
}
\description{
  Method for generating a sampling design for data
  generation following a Binomial-Gaussian model.
}
\details{
  Generates a sampling design for the heterogeneity 'h',
  balancing factors 'a1', ..., 'ak' of group assignments,
  and study sizes 's1', ..., 'sk'.  This design can be used
  for testing methods for inference for the random effects
  meta regression model, since the logarithm of relative
  risks of each study is approximately Gaussian
  distributed. One may use methods that adjust for
  uncertainty in the heteroscedasticity estimates by
  additionally considering the size of the respected
  studies.

  Points in the design are selected via a maxi-min
  hypercube sampling using the 'lhs' package in a
  predefined parameter cube.
}
\examples{
dB <- designB(n=15L, h_bounds=c(0,1), a_bounds=c(-.3,3),
  s_bounds=c(200L,2000L), r=0.03, x=cbind(1,1:5))

if(!all(dim(dB) == c(15,2*dim(cbind(1,1:5))[1]+2))) {
  stop("Wrong dimension")
}
}

