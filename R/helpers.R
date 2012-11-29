#file: ./helpers.R
#author: Thomas Friedrich (info@suud.de)

### The q_delta(tau) function.
###

qfunc <- function (  y # study responses
                   , d # heteroscedasticity
                   , x # design matrix
                   ) {
    ## Returns the q-function. (fastest)
    SVD <- svd(t(x),nv=dim(x)[1])
    K   <- (SVD$v)[,-(1:length(SVD$d))]
    tmp <- y%*%K
    qfuncuno <- function (h) {
        return(as.vector(tmp %*% ginv(t(K) %*% diag(h+d) %*% K) %*% t(tmp)))
    }
    return(function (vec) {return(sapply(vec, qfuncuno))})
}

### The p_delta(eta) function.
###

pfunc <- function ( y   # study responses
                  , d   # heteroscedasticity
                  , x   # design matrix
                  ) {
    ## Returns the p-function.
    qfun  <- qfunc(y,d,x)
    qfu0  <- qfun(0)
    tmpf  <- function(tmp, eta) {return(qfun(tmp) - eta)}
    findk <- 0
    pfun  <- function(eta) {
        findk <- 0
        while (tmpf(exp(findk), eta) > 0) {findk = findk+1}
        return(uniroot(tmpf, c(0,exp(findk)), f.lower=qfu0, eta=eta)$root)}
    pfuncuno <- function(eta) {
        tmp=Inf
        if (eta >= qfu0) {tmp=0} else if (eta > 0) {tmp=pfun(eta)}
        return(tmp)
    }
    return(function(etavec) {return(sapply(etavec, pfuncuno))})
}

### Random streams of pivotal quantities
###

pivotalStream_1 <- function (  n # number of draws for each simulation
                             , y # study response vector
                             , d # heteroscedasticity vector
                             , x # design matrix
                             ) {
    ## Simulate random draws of generalised pivotal quantities.
    k <- dim(x)[1]
    p <- dim(x)[2]
    ## Generating the H-pivot
    pivh <- data.frame(pivh=pfunc(y,d,x)(rchisq(n,df=k-p)))
    ## Generating the N-pivot
    ns  <- matrix(  rnorm(n*p)
                  , ncol=p
                  , dimnames=list(NULL,paste("n",1:p,sep="")))
    return(cbind(ns,pivh))
}

pivotalStream_2 <- function (  n # number of draws for each simulation
                             , y # study response vector
                             , d # heteroscedasticity vector
                             , s # study sizes
                             , x # design matrix
                             ) {
    ## Simulate random draws of generalised pivotal quantities.
    k <- dim(x)[1]
    p <- dim(x)[2]
    ## Generating the D-pivot
    pivd <- d * (s-1) / matrix(  rchisq(n*k, df=s)
                               , ncol=k, byrow=TRUE
                               , dimnames=list(NULL,paste("pivd",1:k,sep="")))
    ## Generating the H-pivot
    makepivh <- function (l) {
        pc <- l[1]
        pd <- l[-1]
        return(pfunc(y,pd,x)(pc))
    }
    pivh <- data.frame(pivh=apply(cbind(rchisq(n,df=k-p), pivd), 1, makepivh))
    ## Generating the N-pivot
    ns  <- matrix(  rnorm(n*p)
                  , ncol=p
                  , dimnames=list(NULL,paste("n",1:p,sep="")))
    return(cbind(ns,pivd,pivh))
}

### Calculate a point estimate from point estimates of variance components
###

calculatePointEstimate <- function (  y # k-vector of responses
                                    , d # k-vector of heteroscedasticities
                                    , h # scalar of heterogeneity
                                    , x # design k-p-matrix
                                    ) {
    ## Calculate a pivotal quantity for the regression coefficents
    ## V = x' * diag(1/(h+d)) * x
    ## B = V^(-1) * x' * diag(1/(h+d))
    O <- diag(1/(h+d))
    V <- t(x) %*% O %*% x
    return(as.vector(solve(V, t(x) %*% O) %*% y))
}

### Calculate pivotal quantities for the regression coefficents
###

makeMultPivot <- function (  y # k-vector of responses
                           , d # k-vector of heteroscedasticities
                           , h # scalar of heterogeneity
                           , n # p-vector of some p-variate Gaussian draw
                           , x # design k-p-matrix
                           ) {
    ## Calculate a pivotal quantity for the regression coefficents
    ## V = x' * diag(1/(h+d)) * x
    ## B = V^(-1) * x' * diag(1/(h+d))
    O <- diag(1/(h+d))
    V <- t(x) %*% O %*% x
    B <- solve(V, t(x) %*% O)
    ## MULTIVARIATE FORMULA
    ## R = (By) - V^(-1/2) * n
    V_eig <- eigen(V)
    evecs <- t(V_eig$vectors)
    R <- as.vector(solve(evecs, evecs%*%B%*%y - diag(sqrt(V_eig$values))%*%evecs%*%n))
    return(R)
}

makeUnivPivot <- function (  y # k-vector of responses
                           , d # k-vector of heteroscedasticities
                           , h # scalar of heterogeneity
                           , n # p-vector of some p-variate Gaussian draw
                           , x # design k-p-matrix
                           ) {
    ## Calculate a pivotal quantity for the regression coefficents
    ## V = x' * diag(1/(h+d)) * x
    ## B = V^(-1) * x' * diag(1/(h+d))
    O <- diag(1/(h+d))
    V <- t(x) %*% O %*% x
    B <- solve(V, t(x) %*% O)
    ## UNIVARIATE FORMULA
    ## L = (By) - sqrt(diag(V^(-1))) * n
    V_simple <- diag(sqrt(diag(ginv(V))))
    L <- as.vector(B%*%y - V_simple %*% n)
    return(L)
}

rRegressionCoefficents_1 <- function (  y         # k-vector of responses
                                      , d         # k-vector of heteroscedasticities
                                      , x         # design k-p-matrix
                                      , rstream   # stream of pivots for the variance components
                                      , rfunction # function to produce pivots for the regCoeffs
                                      ) {
    ## Turns a random stream of pivotal quantities into
    ## a random draws of a pivotal quantity of the
    ## regression coefficents
    pivn <- grep("^n[0-9]", names(rstream))
    pivh <- grep("^pivh", names(rstream))
    makePivots <- function (r) {
        regCoeffPivot <- rfunction(y=y, d=d, h=r[pivh], n=r[pivn], x)
        data.frame(  parameter=1:dim(x)[2]
                   , value=regCoeffPivot)
    }
    return(adply(as.matrix(rstream), 1, makePivots))
}

rRegressionCoefficents_2 <- function (  y    # k-vector of responses
                                      , x    # design k-p-matrix
                                      , rstream   # stream of pivots for the variance components
                                      , rfunction # function to produce pivots for the regCoeffs
                                      ) {
    ## Turns a random stream of pivotal quantities into
    ## a random draws of a pivotal quantity of the
    ## regression coefficents
    pivn <- grep("^n[0-9]", names(rstream))
    pivh <- grep("^pivh", names(rstream))
    pivd <- grep("^pivd[0-9]", names(rstream))
    makePivots <- function (r) {
        regCoeffPivot <- rfunction(y=y, d=r[pivd], h=r[pivh], n=r[pivn], x)
        data.frame(  parameter=1:dim(x)[2]
                   , value=regCoeffPivot)
    }
    return(adply(as.matrix(rstream), 1, makePivots))
}

### Functions that produce generalised confidence intervals
### and point estimates from independend draws of
### generalised pivotal quantities of the regression
### coefficents.
###

makePointEstimates <- function (  rD # pivotal quantities for the regression coefficent vector
                                ) {
    ## rD should be a data frame of values with variabels:
    ## parameter, type and value.  Returns point estimates for
    ## each parameter and type based on the median and the
    ## mean of the pivotal draws
    return(ddply(  rD, "parameter"
                 , function (r) {
                     return(data.frame(  type=factor(c("Generalised Median", "Generalised Expected"))
                                       , value=c(median(r$value), mean(r$value))))
                 }))
}

makeInfByVec <- function ( r    # pivotal quantities for the regression coefficent vector
                         , sgnf # vector of significance levels
                         ) {
    quantiles <- unname(quantile(r, c(sgnf/2, 1-sgnf/2)))
    return(data.frame(  confidence=1-sgnf
                      , lower=quantiles[1:length(sgnf)]
                      , upper=quantiles[-(1:length(sgnf))]))
}

makeInfByDat <- function (  rD # pivotal quantities for the regression coefficent vector
                          , sgnf
                          ) {
    ## rD should be a data frame of values with covariabels
    ## parameter, type and value.  Returns confidence intervals for
    ## each parameter and type.
    return(ddply(  rD, "parameter"
                 , function (r) {return(makeInfByVec(r=r$value, sgnf=sgnf))}))
}
