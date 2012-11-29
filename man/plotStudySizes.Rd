\name{plotStudySizes}
\alias{plotStudySizes}
\title{Example: Plotting study sizes}
\usage{
  plotStudySizes(dat)
}
\arguments{
  \item{dat}{data frame of study responses of binomial
  type.}
}
\value{
  An object created by ggplot2.
}
\description{
  Example: Plotting study sizes
}
\examples{
bcg <- bcgVaccineData()
plotStudySizes(bcg)
}

