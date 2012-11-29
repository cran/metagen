\name{hEstimates}
\alias{hEstimates}
\title{Point estimates: For the heterogeneity parameter}
\usage{
  hEstimates(y, d, x)
}
\arguments{
  \item{y}{study responses}

  \item{d}{heteroscedasticity}

  \item{x}{design matrix}
}
\value{
  A data frame containing point estimates.  Variables are
  'type' and 'h'.
}
\description{
  Returns a list of tau estimates based on different
  approximative methods. Different point estimates for the
  heterogeneity parameter are calculated: HD (Hedges), SL
  (DerSimonian-Laird), SJ (Sidik-Jonkman), MP
  (Mandel-Paule), ML (maximum likelihood), REML (restricted
  maximum-likelihood). Since any of these methods may fail
  to converge, there result may be 'NA' in this case.
}
\examples{
bcg   <- bcgVaccineData()
bcg_y <- bcg$logrisk
bcg_d <- bcg$sdiv
bcg_x <- cbind(1,bcg$x)
hEstimates(y=bcg_y, d=bcg_d, x=bcg_x)

# The implementation can also handle the case in which
# a meta regression is planed with no intercept and only a
# single covariate (i.e. dim(x) = 1).  In this case,
# the design matrix can simply be provided by a vector.
# (This makes no sense in this example and shall only prove
# feasibility)
hEstimates(y=bcg_y, d=bcg_d, x=bcg$x)

# When performing a meta analysis, provide the function
# with a vector of 1s.
hEstimates(y=bcg_y, d=bcg_d, x=rep(1, length(bcg_y)))
}

