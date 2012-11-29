\name{pivotalStream}
\alias{pivotalStream}
\title{Steams of pivotal quantities of the regression coefficient}
\usage{
  pivotalStream(n, y, d, x, s = NULL,
    method = list("univariate", "multivariate"), adjusted)
}
\arguments{
  \item{n}{length of stream.}

  \item{y}{k-vector of responses.}

  \item{d}{k-vector of heteroscedasticity.}

  \item{x}{design (k,p)-matrix.}

  \item{s}{k-vector of study responses.  No need to provide
  this, when adjusted=FALSE.  Default is NULL.}

  \item{method}{A list.  Used to choose the methods for
  calculating the pivotal quantities of the regression
  coefficients.  Default is 'method=list("univariate",
  "multivariate")'.}

  \item{adjusted}{TRUE or FALSE.  Default is FALSE.}
}
\value{
  If method=="univariate" or method=="multivariate", then
  the return is a (p+1)-n-matrix.  The first row contains
  pivotal quantities of the heterogeneity, the rest of the
  rows pivotal quantities of the regression coefficients.
  Each column is an independent draw.

  If 'method==list("univariate", "multivariate")', then the
  return is a (2p+1)-n-matrix.  Of each column, the first
  element is a pivotal for the heterogeneity, the next 'p'
  elements is a pivotal vector for the regression
  coefficients based on "univariate", the last 'p' elements
  are a pivotal vector for the regression coefficients
  based on "multivariate"
}
\description{
  Algorithm for generating a steam of generalised pivotal
  quantities for the regression coefficients.  If
  adjusted=FALSE, then not adjusting for uncertainty in the
  heteroscedasticity estimates d. If adjusted=TRUE, then
  adjusted.  In this case, 's' need to be provided.
}

