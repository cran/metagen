\name{plotDensityH}
\alias{plotDensityH}
\title{Pivotal distributions: Plot pivotal distribution of heterogeneity}
\usage{
  plotDensityH(pivh)
}
\arguments{
  \item{pivh}{pivotal stream with or without adjustment of
  independent draws of a pivotal quantity of the
  heterogeneity.}
}
\description{
  Pivotal distributions: Plot pivotal distribution of
  heterogeneity
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_s <- bcg$size
bcg_x <- cbind(1,bcg$x)

set.seed(865287113)
pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
  adjusted=FALSE)
set.seed(865287113)
pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
  adjusted=TRUE)

pivh <- joinPivotalHeterogeneity(pivUn, pivAd)
plotDensityH(pivh)
}

