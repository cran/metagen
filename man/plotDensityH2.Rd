\name{plotDensityH2}
\alias{plotDensityH2}
\title{Pivotal distributions: Plot pivot density of the heterogeneity}
\usage{
  plotDensityH2(pivh)
}
\arguments{
  \item{pivh}{pivotal stream with or without adjustment of
  independent draws of a pivotal quantity of the
  heterogeneity.}
}
\description{
  Pivotal distributions: Plot pivot density of the
  heterogeneity
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)

set.seed(865287113)
pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
  adjusted=FALSE)
pivh  <- joinPivotalHeterogeneity(pivUn)
plotDensityH2(pivh)
}

