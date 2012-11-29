\name{plotStudyQfuncPfunc}
\alias{plotStudyQfuncPfunc}
\title{Example: Plotting the q- and p-function from the dissertation}
\usage{
  plotStudyQfuncPfunc(y, d, x, n)
}
\arguments{
  \item{y}{a vector of responses.}

  \item{d}{a vector of heteroscedasticity.}

  \item{x}{a design matrix.}

  \item{n}{number of points to interpolate along.}
}
\value{
  A list of objects created by ggplot2.
}
\description{
  Example: Plotting the q- and p-function from the
  dissertation
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_s <- bcg$size
bcg_x <- cbind(1,bcg$x)
p <- plotStudyQfuncPfunc(y=bcg_y, d=bcg_d, x=bcg_x, n=500)
p[1] # plot of the q-function
p[2] # plot of the p-funciton
}

