\name{makeConfInts}
\alias{makeConfInts}
\title{Interval estimates: Generic function}
\usage{
  makeConfInts(sgn, pst, fct, crt, name)
}
\arguments{
  \item{sgn}{one significance level.}

  \item{pst}{point estimate.}

  \item{fct}{standard error.}

  \item{crt}{function for critical value computation.}

  \item{name}{string: name of the method.}
}
\description{
  Generic function to produce interval estimates of
  univariate parameters based on first order limit theory.
}
\details{
  Function for symmetric confidence intervals based on
  standard deviations, point estimates, and quantile
  functions.
}

