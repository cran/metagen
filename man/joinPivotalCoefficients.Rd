\name{joinPivotalCoefficients}
\alias{joinPivotalCoefficients}
\title{Pivotal distributions: Extract pivots for regression coefficients}
\usage{
  joinPivotalCoefficients(p0, p1)
}
\arguments{
  \item{p0}{pivotal stream without adjustment}

  \item{p1}{pivatal stream with adjustment}
}
\description{
  Pivotal distributions: Extract pivots for regression
  coefficients
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

pivr <- joinPivotalCoefficients(pivUn, pivAd)
}

