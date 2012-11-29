# Copyright (C) 2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

#' Pivotal distributions: Extract pivots for heterogeneity
#'
#' @param p0 pivotal stream without adjustment
#' @param p1 pivatal stream with adjustment
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivh <- joinPivotalHeterogeneity(pivUn, pivAd)
#' @export
joinPivotalHeterogeneity <- function (p0=NULL, p1=NULL) {
    if (!is.null(p0)) {
        p0 <- data.frame(  heterogeneity=p0[1,]
                         , type=factor(  1
                                       , levels=1:2
                                       , labels=c("unadjusted",
                                                  "adjusted")
                                       , ordered=TRUE))
    }

    if (!is.null(p1)) {
        p1 <- data.frame(  heterogeneity=p1[1,]
                         , type=factor(  2
                                       , levels=1:2
                                       , labels=c("unadjusted",
                                                  "adjusted")
                                       , ordered=TRUE))
    }
    return(rbind(p0, p1))
}

getPivotalData <- function(p, vec, name) {
    p <- t(p[vec,])
    colnames(p) <- c("intercept", "slope")
    p <- as.data.frame(p)
    p$via <- name
    return(p)
}

makePivotalData <- function (p) {
    pivRL <- rbind(  getPivotalData(p, 2:3, 1)
                   , getPivotalData(p, 4:5, 2))
    pivRL$via  <- with(  pivRL
                       , factor(  via
                                , levels=2:1
                                , labels=c("formula R", "formula L")
                                , ordered = TRUE))
    return(pivRL)
}

#' Pivotal distributions: Extract pivots for regression coefficients
#'
#' @param p0 pivotal stream without adjustment
#' @param p1 pivatal stream with adjustment
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#'   pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#' adjusted=TRUE)
#'
#' pivr <- joinPivotalCoefficients(pivUn, pivAd)
#' @export
joinPivotalCoefficients <- function (p0, p1) {
    pivUn <- makePivotalData(p0)
    pivAd <- makePivotalData(p1)

    pivUn$type <- 1
    pivAd$type <- 2
    pivDat <- rbind(pivUn, pivAd)
    pivDat$type <- with(pivDat, factor(  type
                                       , levels=1:2
                                       , labels=c("unadjusted",
                                                  "adjusted")
                                       , ordered=TRUE))
    return(pivDat)
}

#' Pivotal distributions: Plot pivotal distribution of heterogeneity
#'
#' @param pivh pivotal stream with or without adjustment of
#' independent draws of a pivotal quantity of the heterogeneity.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivh <- joinPivotalHeterogeneity(pivUn, pivAd)
#' plotDensityH(pivh)
#' @export
plotDensityH <- function (pivh) {
    return( ggplot(pivh, aes(x=heterogeneity))
           + geom_density(aes(fill=type), alpha=.42)
           + xlab(expression(tau))
           + scale_fill_brewer(  expression(atop(  "Uncertainty in the"
                                                 , "heteroscedasticity"
                                                 , "has been"))
           , palette="Spectral")
           + ggtitle(expression("Density estimates of generalised pivots
                                for the heterogeneity"))
           )
}

#' Pivotal distributions: Plot pivot density of the heterogeneity
#'
#' @param pivh pivotal stream with or without adjustment of
#' independent draws of a pivotal quantity of the heterogeneity.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' pivh  <- joinPivotalHeterogeneity(pivUn)
#' plotDensityH2(pivh)
#' @export
plotDensityH2 <- function (pivh) {
    return(ggplot(pivh, aes(x=heterogeneity))
           + geom_histogram(aes(y=..density.., fill=..count..),
                            alpha=.65)
           + geom_density()
           + scale_fill_gradient(expression(atop(  "Counts within the"
                                                 , "respected bin"))
                                 , low="darkblue", high="darkgreen"
                                 )
           + xlab(expression(tau))
           + facet_wrap(~type, ncol=1)
           + ggtitle(expression(strwrap("Density estimates of a
                                        generalised pivot for the
                                        heterogeneity")))
           )
}

#' Pivotal distributions: Plot pivotal distribution of regression
#' coefficients
#'
#' @param pivr data frame of independent draws from of pivots.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivr <- joinPivotalCoefficients(pivUn, pivAd)
#' plotDensityIntercept(pivr)
#' @export
plotDensityIntercept <- function (pivr) {
    return( ggplot(pivr, aes(x=intercept))
           + geom_density(aes(fill=via), alpha=.42)
           + xlab(expression("intercept"))
           + scale_fill_brewer(palette="Spectral")
           + facet_grid(~type, scales="free")
           + ggtitle(expression(strwrap("Density estimates of
                                        generalised pivots for the
                                        regression intercept")))
           )
}

#' Pivotal distributions: Plot pivotal distribution of regression
#' coefficients
#'
#' @param pivr data frame of independent draws from of pivots.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivr <- joinPivotalCoefficients(pivUn, pivAd)
#' plotDensitySlope(pivr)
#' @export
plotDensitySlope <- function (pivr) {
    return( ggplot(pivr, aes(x=slope))
           + geom_density(aes(fill=via), alpha=.42)
           + xlab(expression("slope"))
           + scale_fill_brewer(  expression("using")
                               , palette="Spectral")
           + facet_grid(~type, scales="free")
           + ggtitle(expression(strwrap("Density estimates of
                                        generalised pivots for the
                                        regression slope")))
           )
}

#' Pivotal distributions: Plot pivotal distribution of regression
#' coefficients
#'
#' @param pivr data frame of independent draws from of pivots.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivr <- joinPivotalCoefficients(pivUn, pivAd)
#' plotDensityIntercept2(pivr)
#' @export
plotDensityIntercept2 <- function (pivr) {
    return( ggplot(pivr, aes(x=intercept))
           + geom_density(aes(fill=type), alpha=.42)
           + xlab(expression("slope"))
           + scale_fill_brewer(  expression(atop(  "Uncertainty in the"
                                                 , "heteroscedasticity"
                                                 , "has been"))
                               , palette="Spectral")
           + facet_grid(~via, scales="free")
           + ggtitle(expression(strwrap("Density estimates of
                                        generalised pivots for the
                                        regression slope")))
           )
}

#' Pivotal distributions: Plot pivotal distribution of regression
#' coefficients
#'
#' @param pivr data frame of independent draws from of pivots.
#' @examples
#' bcg   <- bcgVaccineData()
#' bcg_y <- bcg$logrisk
#' bcg_d <- bcg$sdiv
#' bcg_s <- bcg$size
#' bcg_x <- cbind(1,bcg$x)
#'
#' set.seed(865287113)
#' pivUn <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x,
#'   adjusted=FALSE)
#' set.seed(865287113)
#' pivAd <- pivotalStream(50, y=bcg_y, d=bcg_d, x=bcg_x, s=bcg_s,
#'   adjusted=TRUE)
#'
#' pivr <- joinPivotalCoefficients(pivUn, pivAd)
#' plotDensitySlope2(pivr)
#' @export
plotDensitySlope2 <- function (pivr) {
    return( ggplot(pivr, aes(x=slope))
           + geom_density(aes(fill=type), alpha=.42)
           + xlab(expression("slope"))
           + scale_fill_brewer(  expression(atop(  "Uncertainty in the"
                                                 , "heteroscedasticity"
                                                 , "has been"))
                               , palette="Spectral")
           + facet_grid(~via, scales="free")
           + ggtitle(expression(strwrap("Density estimates of
                                        generalised pivots for the
                                        regression slope")))
           )
}

globalVariables()
