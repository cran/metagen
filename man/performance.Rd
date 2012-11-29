\name{performance}
\alias{performance}
\title{Running a computer experiment}
\usage{
  performance(results, b, h)
}
\arguments{
  \item{results}{Needs to be of the same type as, for
  example, the return value of the computer experiments
  'experimentY', 'experimentD'.}

  \item{b}{true regression coefficients}

  \item{h}{true heterogeneity}
}
\value{
  Data frame containing performance measuments of inference
  methods based on the results of the computer experiment
  given by 'results'.
}
\description{
  Adding performance measures to the results
}
\details{
  Calculating performance measurements from a computer
  experiment for analysis data from a random effects meta
  regression model.
}
\examples{
h_test <- 0.03
x_test <- cbind(1,1:7)
b_test <- c(.5, .25)
sgnf_test <- c(0.025, 0.01)

set.seed(5133568) # for reproducibility
d_test <- rchisq(7, df=0.02)

# In an actual computer experiment, use 'piv_draws=1000' instead!!
eY <- experimentY(n=5, h=h_test, d=d_test, x=x_test, b=b_test,
  sgnf=sgnf_test, piv_draws=50)

performance(results=eY, b=b_test, h=h_test)
}

