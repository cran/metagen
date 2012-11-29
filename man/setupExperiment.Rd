\name{setupExperiment}
\alias{setupExperiment}
\title{Running a computer experiment in batch mode}
\usage{
  setupExperiment(name, seed, n, resolution, bounds, x, b,
    sgnf, piv_draws, ...)
}
\arguments{
  \item{name}{Reference name for the experiment.}

  \item{seed}{Random seed for the experiment.}

  \item{n}{number of simulations to at each parameter
  configuration.}

  \item{resolution}{list of number of parameter
  configurations in each design, e.g.
  resolution=list(h=5L, d=3L)}

  \item{bounds}{list of parameter bounds used for
  experimental design, e.g.  bounds=list(h=c(0,1),
  d=c(0.001, 2), s=c(200L, 2000L)) where - h : bounds of
  the heterogeneity.  - d : bounds of the
  heteroscedasticity.  - a : bounds of the balancing factor
  of group assignments.  - s : bounds of the study sizes.
  - r : fixed risk in the control.}

  \item{x}{design matrix.}

  \item{b}{regression coefficients.}

  \item{sgnf}{levels of significance.}

  \item{piv_draws}{number of pivotal draws.}

  \item{...}{further arguments to makeExperimentRegistry,
  e.g.  file.dir=tempfile().}
}
\value{
  The registry.
}
\description{
  Sets up a computer experiment evaluating the performance
  of different inference methods in the random effects meta
  regression model.
}

