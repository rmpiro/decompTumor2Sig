

###################################
# internal mathematical functions #
###################################

#' computeFrobeniusNorm (internal function)
#'
#' Compute the Frobenius norm of a numeric matrix (if this matrix is the 
#' difference between two matrices, then this corresponds to the Frobenius 
#' distance between these two matrices).
#'
#' @usage computeFrobeniusNorm(A)
#' @param A The numeric matrix.
#' @return The Frobenius norm.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
computeFrobeniusNorm <- function(A) { # compute the Frobenius norm of matrix A
    A <- as.matrix(A) # in case A was a data.frame

    trace <- sum(diag(A %*% t(A)))      # need the conjugate transpose here, 
                                        # but for real-valued matrices this is 
                                        # equivalent to the transpose!

    return(sqrt(trace))
}


#' computeRSS (internal function)
#'
#' Compute the residual sum of squares (RSS), i.e., the sum of squared errors.
#'
#' @usage computeRSS(x, y)
#' @param x The first numeric object (e.g., matrix).
#' @param y The second numeric object (of the same type as \code{x}).
#' @return The RSS.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
computeRSS <- function(x, y) {
    if (length(x) != length(y)) {
        stop("In computeRSS: x and y need to have the same length!")
    }

    return(sum((x-y)^2))
}


#' roundIntegerSum (internal function)
#'
#' update a numeric vector such that it's composed of integers and its sum
#' reaches a desired total. Positive or negative discrepancies are distributed
#' proportionally between the summands.
#'
#' @usage roundIntegerSum(vec, targetSum)
#' @param vec Vector of integers.
#' @param targetSum The target sum to be reached.
#' @return The updated vector with the desired total sum.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
roundIntegerSum <- function(vec, targetSum) {
    stopifnot(is.numeric(vec) || is.integer(vec))
    stopifnot(is.numeric(targetSum) || is.integer(targetSum))

    # first get an approximation
    vec <- round(targetSum * (vec / sum(vec)))
    
    # modification order: highest first, etc.
    modOrder <- order(vec, decreasing=TRUE)

    # See if we still have some discrepancy
    diff <- targetSum - sum(vec)    # pos: need to add; neg: need to subtract

    # The difference should be between -1 and 1, not beyond this range
    if (diff == -1) {
        vec[modOrder[length(vec)]] <- vec[modOrder[length(vec)]] - 1
    } else if (diff == 1) {
        vec[modOrder[1]] <- vec[modOrder[1]] + 1
    } else {
        stopifnot(!diff)
    }

    vec
}
