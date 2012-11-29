#' Generalised Inference in the Random Effects Meta Regression Model
#'
#' Provides methods for making inference in the random effects meta
#' regression model, such as point estimates and confidence intervals
#' for the heterogeneity parameter and the regression coefficients
#' vector.
#'
#' The function metagen assumes the model
#'
#'   Y ~ N_k(xb, [h]_k + [d])
#'
#' where Y is a k-vector of study responses, N_k is a k-variate
#' Gaussian distribution, x is a design k-p-matrix, h is a parameter of
#' heterogeneity, [h]_k is the diagonal matrix with h on the diagonal,
#' d is a k-vector of heteroscedasticity and [d] is the diagonal matrix
#' with (d_1,...,d_k) on the diagonal.
#'
#' The function metagen2 assumes the model
#'
#'   Y   ~ N_k(xb, [h]_k + [D])\cr
#'   D_j = K_j d_j\cr
#'   K_j / (s_j - 1) ~ Chisq(s_j-1)
#'
#' where D=(D_1,...,D_k), n_j is the total number of subjects in study j
#' and Chisq(n_j-1) as a Chi-Square distribution with
#' (n_j - 1) degrees of freedom.
#'
#' @param y [\code{numeric(k)}]\cr
#'   Vector of observed study responses.
#' @param d [\code{numeric(k)}]\cr
#'   Vector of heteroscedasticity.
#' @param x [\code{logical(1)}]\cr
#'   The design matrix.
#' @param sgnf [\code{numeric()}]\cr
#'   Significance level.  This can also be a vector of different significance levels.  Please note that you will need to provide a significance level here, not a confidence level.
#' @param numb [\code{integer(1)}]\cr
#'   Number of draws from the distribution of the generalised pivotal quantity.  Usually, 1500-2000 is a good choice here.
#' @return [\code{list}].
#'   \describe{
#'     \item{confh [\code{data.frame}]}{Confidence intervals for the heterogeneity parameter.}
#'     \item{confr [\code{data.frame}]}{Confidence intervals for the regression coefficients.}
#'     \item{pointh [\code{data.frame}]}{Point estimates of the heterogeneity parameter.  One estimate is based on the median, the other is based on the mean of the distribution of the generalised pivotal quantity of the heterogeneity parameter.}
#'     \item{pointr [\code{data.frame}]}{Point estimates of the regression coefficients.  One estimate is based on the median, the other is based on the mean of the distribution of the generalised pivotal quantity of the regression coefficients.}
#' }
#' @export
#' @examples
#' library(metafor)
#' data(dat.bcg)
#' dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, append=TRUE)
#' res <- metagen(y=dat$yi, d=dat$vi, x=cbind(1, dat$ablat))
#'
#' sizes <- dat$tpos + dat$tneg + dat$cpos + dat$cneg
#' res2 <- metagen2(y=dat$yi, d=dat$vi, s=sizes, x=cbind(1, dat$ablat))
metagen <- function (  y    # k-vector of study responses
                     , d    # k-vector of heteroscedasticity
                     , x    # design k-p-matrix
                     , sgnf=c(0.025, 0.05) # significance levels
                     , numb=1500L # number of draws
                     ) {
    checkArg(y, "numeric", min.len=2L, na.ok=FALSE)
    checkArg(d, "numeric", len=length(y), na.ok=FALSE)
    checkArg(x, "matrix", na.ok=FALSE)
    checkArg(sgnf, "numeric", min.len=1L, lower=0, upper=1, na.ok=FALSE)
    numb = convertInteger(numb)
    checkArg(numb, "integer", len=1L, lower=1, na.ok=FALSE)

    rstream <- pivotalStream_1(n=numb, y=y, d=d, x=x)

    ## Inference on the heterogeneity parameter
    hconf  <- makeInfByVec(rstream$pivh, sgnf)
    hpoint <- data.frame(  type=factor(c("Generalised Median", "Generalised Expected"))
                         , h=c(median(rstream$pivh), mean(rstream$pivh)))
    ## Inference on the regression coefficent parameters
    rPivM   <- rRegressionCoefficents_1(y,d,x,rstream,makeMultPivot)
    rConfM  <- cbind(method="Generalised Multivariate", makeInfByDat(rPivM,sgnf))
    rPointM <- cbind(method="Generalised Multivariate", makePointEstimates(rPivM))
    rPivU   <- rRegressionCoefficents_1(y,d,x,rstream,makeUnivPivot)
    rConfU  <- cbind(method="Generalised Univariate", makeInfByDat(rPivU,sgnf))
    rPointU <- cbind(method="Generalised Univariate", makePointEstimates(rPivU))
    return(list(  confh=hconf
                , confr=rbind(rConfM,rConfU)
                , pointh=hpoint
                , pointr=rbind(rPointM,rPointU)))
}

#' @param s [\code{integer(k)}]\cr
#'   Vector of total study sizes.
#' @export
#' @rdname metagen
metagen2 <- function (  y    # k-vector of study responses
                      , d    # k-vector of heteroscedasticity
                      , s    # k-vector of study sizes
                      , x    # design k-p-matrix
                      , sgnf=c(0.025, 0.05) # significance levels
                      , numb=1500L # number of draws
                      ) {
    checkArg(y, "numeric", min.len=2L, na.ok=FALSE)
    checkArg(d, "numeric", len=length(y), na.ok=FALSE)
    s = convertInteger(s)
    checkArg(s, "integer", len=length(y), na.ok=FALSE)
    checkArg(x, "matrix", na.ok=FALSE)
    checkArg(sgnf, "numeric", min.len=1L, lower=0, upper=1, na.ok=FALSE)
    numb = convertInteger(numb)
    checkArg(numb, "integer", len=1L, lower=1, na.ok=FALSE)

    rstream <- pivotalStream_2(n=numb, y=y, d=d, s=s, x=x)

    ## Inference on the heterogeneity parameter
    hconf  <- makeInfByVec(rstream$pivh, sgnf)
    hpoint <- data.frame(  type=factor(c("Generalised Median", "Generalised Expected"))
                         , h=c(median(rstream$pivh), mean(rstream$pivh)))
    rPivM   <- rRegressionCoefficents_2(y,x,rstream,makeMultPivot)
    rConfM  <- cbind(method="Generalised Multivariate", makeInfByDat(rPivM,sgnf))
    rPointM <- cbind(method="Generalised Multivariate", makePointEstimates(rPivM))
    rPivU   <- rRegressionCoefficents_2(y,x,rstream,makeUnivPivot)
    rConfU  <- cbind(method="Generalised Univariate", makeInfByDat(rPivU,sgnf))
    rPointU <- cbind(method="Generalised Univariate", makePointEstimates(rPivU))
    return(list(  confh=hconf
                , confr=rbind(rConfM,rConfU)
                , pointh=hpoint
                , pointr=rbind(rPointM,rPointU)))
}
