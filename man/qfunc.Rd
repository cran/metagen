\name{qfunc}
\alias{qfunc}
\title{The q_delta(tau) function.}
\usage{
  qfunc(y, d, x)
}
\arguments{
  \item{y}{study responses.}

  \item{d}{heteroscedasticity.}

  \item{x}{design matrix.}
}
\value{
  A vector valued function.
}
\description{
  Returns the q-function, in the form of the so-far fastest
  version tested.
}
\details{
  FIXME: Performance may be improved further!
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)
qfunc(y=bcg_y, d=bcg_d, x=bcg_x)
}

