\name{plotStudyForest}
\alias{plotStudyForest}
\title{Example: Plotting a forest plot of a data frame}
\usage{
  plotStudyForest(dat)
}
\arguments{
  \item{dat}{data frame of study responses of binomial
  type.}
}
\value{
  An object created by ggplot2.
}
\description{
  Example: Plotting a forest plot of a data frame
}
\examples{
bcg <- bcgVaccineData()
plotStudyForest(bcg)
}

