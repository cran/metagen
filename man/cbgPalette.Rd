\docType{data}
\name{cbgPalette}
\alias{cbgPalette}
\title{Colour palettes for colour blind people}
\format{chr [1:8] "#999999" "#E69F00" "#56B4E9" "#009E73" "#F0E442" ...}
\usage{
  cbgPalette
}
\description{
  The palette with grey.
}
\details{
  This palette is directly taken from

  http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

  Hence, I don't take any credit for this.
}
\examples{
scale_fill_discrete <- function(...) scale_fill_manual(...,
  values=cbgPalette)
scale_colour_discrete <- function(...) scale_fill_manual(...,
  values=cbgPalette)
}
\keyword{datasets}

