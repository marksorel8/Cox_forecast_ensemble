residuals.dlmFiltered <- function(object, ...,

                                  type=c("standardized", "raw"), sd=TRUE) {

    if (is.null(object$y))

        stop("\'object\' argument has no \'y\' component")

    type <- match.arg(type)

    if (is.null(object$mod$JFF)) tvFF <- FALSE else tvFF <- TRUE

    if (is.null(object$mod$JV)) tvV <- FALSE else tvV <- TRUE

    FF <- object$mod$FF

    if (!( tvFF || tvV )) { ## constant model

        f <- as.matrix(object$a) %*% t(FF)

        res <- drop(object$y - f) # one-step forecasting errors

        if (sd || (type == "standardized")) {

            V <- object$mod$V

            SD <- drop(t(sqrt(sapply(seq(along=object$U.R),

                                     function(i)

                                     diag(crossprod(object$D.R[i,] *

                                                    t(FF%*%object$U.R[[i]])) + V)))))

        }

    } else

    if ( !tvFF ) { ## only V time-varying

        f <- as.matrix(object$a) %*% t(FF)

        res <- drop(object$y - f) # one-step forecasting errors

        if (sd || (type == "standardized")) {

            nz <- object$mod$JV != 0

            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],

                        object$mod$JV[nz])

            V <- object$mod$V

            getSD <- function(i) {

                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]

                diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V)

            }

            SD <- drop(t(sqrt(sapply(seq(along=object$U.R), getSD))))

        }

    } else

    if ( !tvV ) { ## only FF time-varying

        if (!(sd || (type == "standardized"))) {

            nz <- object$mod$JFF != 0

            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],

                         object$mod$JFF[nz])

            getFore <- function(i) {

                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]

                FF %*% as.matrix(object$a)[i,]

            }

            f <- drop(t(sapply(seq(along=object$U.R), getFore)))

            res <- drop(object$y - f) # one-step forecasting errors

        } else {

            nz <- object$mod$JFF != 0

            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],

                         object$mod$JFF[nz])

            V <- object$mod$V

            getBoth <- function(i) {

                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]

                c(FF %*% as.matrix(object$a)[i,],

                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))

            }

            tmp <- t(sapply(seq(along=object$U.R), getBoth))

            m <- ncol(tmp) / 2

            res <- drop(object$y - tmp[,1:m])

            SD <- drop(sqrt(tmp[,-(1:m)]))

        }

    } else { ## both FF and V time-varying

        if (!(sd || (type == "standardized"))) {

            nz <- object$mod$JFF != 0

            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],

                         object$mod$JFF[nz])

            getFore <- function(i) {

                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]

                FF %*% as.matrix(object$a)[i,]

            }

            f <- drop(t(sapply(seq(along=object$U.R), getFore)))

            res <- drop(object$y - f) # one-step forecasting errors

        } else {

            nz <- object$mod$JFF != 0

            JFF <- cbind(row(object$mod$JFF)[nz], col(object$mod$JFF)[nz],

                         object$mod$JFF[nz])

            nz <- object$mod$JV != 0

            JV <- cbind(row(object$mod$JV)[nz], col(object$mod$JV)[nz],

                        object$mod$JV[nz])

            V <- object$mod$V

            getBoth <- function(i) {

                FF[JFF[,-3,drop=FALSE]] <- object$mod$X[i,JFF[,3]]

                V[JV[,-3,drop=FALSE]] <- object$mod$X[i,JV[,3]]

                c(FF %*% as.matrix(object$a)[i,],

                  diag(crossprod(object$D.R[i,] * t(FF%*%object$U.R[[i]])) + V))

            }

            tmp <- t(sapply(seq(along=object$U.R), getBoth))

            m <- ncol(tmp) / 2

            res <- drop(object$y - tmp[,1:m])

            SD <- drop(sqrt(tmp[,-(1:m)]))

        }

    }

 

    if ( type == "standardized" )

        res <- res / SD

    if (sd) {

        if (is.ts(res)) attributes(SD) <- attributes(res) # makes a time series of SD

        return(list(res=res, sd=SD))

    } else

    return(res)

}
