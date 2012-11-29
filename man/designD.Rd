\name{designD}
\alias{designD}
\title{Design: Gaussian responses (unknown heteroscedasticity)}
\usage{
  designD(n, h_bounds, d_bounds, s_bounds, x)
}
\arguments{
  \item{n}{resolution of the heterogeneity and
  heteroscedasticity parameters, i.e., the number of of
  different (heterogeneity, heteroscedasticity, sizes)
  tuple in the design.}

  \item{h_bounds}{bounds of the heterogeneity.}

  \item{d_bounds}{bounds of the heteroscedasticity.}

  \item{s_bounds}{bounds of the study sizes.}

  \item{x}{design matrix.}
}
\value{
  Function returns a data frame.  Each line of this data
  frame can be an input to the function 'rD' which is used
  to sample data from such a design.
}
\description{
  Method for generating a sampling design for data
  generation following a random effects meta regression
  model with unknown heteroscedasticity.
}
\details{
  Generates a sampling design for the heterogeneity 'h',
  heteroscedasticity 'd1', ..., 'dk', and study sizes 's1',
  ..., 'sk'. This design can be used for testing methods
  that adjust for uncertainty in the heteroscedasticity
  estimates by additionally considering the size of the
  respected studies.

  Points in the design are selected via a maxi-min
  hypercube sampling using the 'lhs' package in a
  predefined parameter cube.
}
\examples{
dD <- designD(n=15L, h_bounds=c(0,1), d_bounds=c(0.01,2),
  s_bounds=c(200L,2000L), x=cbind(1,1:7))

if(!all(dim(dD) == c(15,2*dim(cbind(1,1:7))[1]+1))) {
  stop("Wrong dimension")
}
}

