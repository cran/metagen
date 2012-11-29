\name{joinPivotalHeterogeneity}
\alias{joinPivotalHeterogeneity}
\title{Pivotal distributions: Extract pivots for heterogeneity}
\usage{
  joinPivotalHeterogeneity(p0 = NULL, p1 = NULL)
}
\arguments{
  \item{p0}{pivotal stream without adjustment.}

  \item{p1}{pivatal stream with adjustment.}
}
\description{
  Pivotal distributions: Extract pivots for heterogeneity
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
}

