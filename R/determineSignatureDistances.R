#' Determine differences between a given signature and a set of
#' target signatures.
#'
#' `determineSignatureDistances()` determines all distances (i.e., differences)
#' between a given signature (of type Alexandrov or Shiraishi) and a set
#' of target signatures (of the same type). This can help to compare signatures
#' that have been determined in different ways or from different datasets.
#' Different distance measures can be used (see details below).
#'
#' Distances that can be used are:
#' 
#' \tabular{ll}{
#' \code{"frobenius"} \tab Forbenius distance between real-valued matrices\cr
#'                    \tab (or Shiraishi signatures) \code{A} and \code{B}:\cr
#'                    \tab \code{F = sqrt(trace( (A-B) \%*\% t(A-B) ))} \cr
#' \code{"rss"}       \tab Residual sum of squares (i.e., squared error):\cr
#'                    \tab \code{rss = sum((A-B)^2)} \cr
#' \code{"euclidean"} \tab (see \code{?dist} for details)\cr
#' \code{"maximum"}   \tab (see \code{?dist} for details)\cr
#' \code{"manhattan"} \tab (see \code{?dist} for details)\cr
#' \code{"canberra"}  \tab (see \code{?dist} for details)\cr
#' \code{"binary"}    \tab (see \code{?dist} for details)\cr
#' \code{"minkowski"} \tab (see \code{?dist} for details)\cr
#' }
#' 
#' @usage determineSignatureDistances(fromSignature, toSignatures,
#'                                    method="euclidean")
#' @param fromSignature (Mandatory) A single signature of the Alexandrov
#' (vector) or Shiraishi type (data frame or matrix). 
#' @param toSignatures (Mandatory) The list of target signatures for which to
#' compute the distances to \code{fromSignature}. These target signatures
#' must be of the same type and format as \code{fromSignature}. 
#' @param method (Optional) The distance measure to be used. This can be one
#' of the following: \code{"frobenius"} for Frobenius distance between
#' matrices (only for Shiraishi signatures); \code{"rss"} for the residual
#' sum of squares (squared error); or any distance measure available for
#' the function \code{dist} of the \code{stats} package.
#' Default: \code{"euclidean"}.
#' @return A signature-named ve ctor containing all distances. This vector has
#' the same order as the target signature list, so it is not sorted according
#' to distance.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{mapSignatureSets}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' ### convert them to Shiraishi signatures
#' signAlex2Shi <- convertAlexandrov2Shiraishi(signAlexandrov)
#' 
#' ### define an arbitrary signature just for testing
#' ### (similar to signature 1)
#' testSig <- matrix(c(0.1,   0, 0.7, 0.1, 0.1,   0,
#'                     0.3, 0.2, 0.3, 0.2,   0,   0,
#'                     0.2, 0.1, 0.5, 0.2,   0,   0), nrow=3, byrow=TRUE) 
#' 
#' ### compute distances of the test signature to the converted
#' ### Alexandrov signatures from COSMIC
#' determineSignatureDistances(testSig, signAlex2Shi, method="frobenius")
#' 
#' @importFrom stats dist
#' @export determineSignatureDistances
determineSignatureDistances <- function(fromSignature, toSignatures,
                                        method="euclidean") {

    if (!isSignatureSet(toSignatures)) {
        stop("Parameter 'toSignatures' must be a list of signature objects!")
    }
    
    # is the toSignatures are unnamed, name them by enumerating them
    if(is.null(names(toSignatures))) {
        names(toSignatures) <- paste0("sign_",seq_along(toSignatures))
    }

    if (!is.probability.object(fromSignature)) {
        
        stop(paste("'fromSignature' must be a single Alexandrov signature",
                   "(data.frame or matrix) or Shiraishi signature (vector)!"))
    }

    if (!sameSignatureFormat(list(fromSignature), toSignatures)) {
        stop("fromSignature and toSignatures must be of the same format.")
    }

    
    # allow Frobenius distance only for data.frame or matrix
    if (method == "frobenius" &
        !(isShiraishiSet(list(fromSignature))) ) {

        stop(paste("Frobenius distance can be used only for Shiraishi",
                   "signatures (matrix or data.frame)!"))
    }
    
    distances <- vapply(toSignatures, function(s) {
        # first, check we have the same format as fromSignature
        if (as.vector(length(fromSignature)) != as.vector(length(s))) {
            stop("Formats of fromSignature and toSignatures must match!")
        }

        # determine distance
        if (method == "frobenius") {
            # Frobenius distance? use own implementation
            d <- computeFrobeniusNorm(as.matrix(fromSignature)-as.matrix(s))
        } else if (method == "rss") {
            # Residual sum of squares (squared error)
            d <- computeRSS(as.vector(fromSignature),as.vector(s))
        } else {
            # other distances? use "dist" from "stats"
            d <- dist(rbind(as.vector(fromSignature),as.vector(s)),
                      method=method)
        }

        d # return the computed distance
    }, FUN.VALUE=numeric(1) )
    
    #names(distances) <- names(toSignatures)
    #sort(distances)

    return(distances)
}

