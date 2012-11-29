\name{bcgVaccineData}
\alias{bcgVaccineData}
\title{Example: Setting up the BCG-data set}
\usage{
  bcgVaccineData(sgnf = 0.05)
}
\arguments{
  \item{sgnf}{significance level of the confidence
  intervals for the relative risks.}
}
\value{
  Returns a data set of 13 clinical trials which evaluated
  the efficacy of the BCG vaccine.  The data set is an
  exact copy of the data set found in the dat.bcg data
  frame provided by the metafor package.
}
\description{
  Exemplary data set of 14 clinical trials evaluating BCG
  vaccine efficacy.
}
\details{
  Reads in the BCG vaccine efficacy data from the metafor
  package and adds some statistics to the data such as the
  log-relative risk, study size, measurements of balance,
  confidence intervals of the responses, and the like.
}
\examples{
bcgVaccineData()
}

