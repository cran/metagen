\name{designY}
\alias{designY}
\title{Design: Gaussian responses (known heteroscedasticity)}
\usage{
  designY(n, h_bounds, d_bounds, x)
}
\arguments{
  \item{n}{resolution of the heterogeneity and
  heteroscedasticity parameters, i.e. the number of of
  different (heterogeneity, heteroscedasticity) pairs in
  the design.}

  \item{h_bounds}{bounds of the heterogeneity.}

  \item{d_bounds}{bounds of the heteroscedasticity.}

  \item{x}{design matrix.}
}
\value{
  Function returns a data frame.  Each line of this data
  frame can be an input to the function 'rY' which is used
  to sample data from such a design.
}
\description{
  Method for generating a sampling design for data
  generation following a random effects meta regression
  model with known heteroscedasticity.
}
\details{
  Generates a sampling design for the heterogeneity 'h' and
  a heteroscedasticity 'd1', ..., 'dk'.

  Points in the design are selected via a maxi-min
  hypercube sampling using the 'lhs' package in a
  predefined parameter cube.
}
\examples{
dY <- designY(n=15L, h_bounds=c(0,1), d_bounds=c(0.01,2),
x=cbind(1,1:7))

if(!all(dim(dY) == c(15,dim(cbind(1,1:7))[1]+1))) {
  stop("Wrong dimension")
}
}

