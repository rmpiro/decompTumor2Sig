#' Map one signature set to another.
#'
#' `mapSignatureSets()` determines a mapping from one set of signatures to
#' another. Both Alexandrov and Shiraishi signatures can be handled, but both
#' sets must be of the same type. The mapping can either be a unique
#' (one-to-one) mapping or identify best matches while allowing multiple
#' signatures to be mapped to the same target signature if it is the best
#' match for more than one signature.
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
#' @usage mapSignatureSets(fromSignatures, toSignatures, method="euclidean",
#' unique=FALSE)
#' @param fromSignatures (Mandatory) A set (list) of signatures of the
#' Alexandrov (vector) or Shiraishi type (data frame or matrix), that has
#' to be mapped to the signatures of a second set (\code{toSignatures}).
#' @param toSignatures (Mandatory) The set (list) of signatures to which the
#' set of \code{fromSignatures} has to be mapped.
#' @param method (Optional) The distance measure to be used. This can be one
#' of the following: \code{"frobenius"} for Frobenius distance between matrices
#' (only for Shiraishi signatures); \code{"rss"} for the residual sum of squares
#' (squared error); or any distance measure available for the function
#' \code{dist()} of the \code{stats} package. Default: \code{"euclidean"}.
#' @param unique (Optional) If set to \code{FALSE} (default), then for each
#' signature of \code{fromSignatures} the best match (minimum distance) from
#' \code{toSignatures} is selected. The selected signatures need not be unique,
#' i.e., one signature of \code{toSignatures} may be the best match for
#' multiple signatures of \code{fromSignatures}. If set to \code{TRUE}, i.e.,
#' if a unique (one-to-one) mapping is required, an iterative approach is
#' performed: in each step, the best matching pair from \code{fromSignatures}
#' and \code{toSignatures} is mapped and then removed from the list of
#' signatures that remain to be mapped, such that they cannot be selected
#' again.
#' @return A vector having as elements the mapped signatures of
#' \code{toSignatures}, and as names the signatures of \code{fromSignatures}
#' with which they have been associated.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{determineSignatureDistances}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' ### convert them to Shiraishi signatures
#' signAlex2Shi <- convertAlexandrov2Shiraishi(signAlexandrov)
#' 
#' ### define a small set of arbitrary signatures just for testing
#' ### (similar to signatures 1, 5 and 13, respectively)
#' test1 <- matrix(c( 0.1,  0,    0.7,  0.1,  0.1,  0,
#'                    0.3,  0.2,  0.3,  0.2,  0,    0,
#'                    0.2,  0.1,  0.5,  0.2,  0,    0   ), nrow=3, byrow=TRUE)
#' 
#' test2 <- matrix(c( 0.1,  0.1,  0.3,  0.1,  0.3,  0.1,
#'                    0.3,  0.25, 0.2,  0.25, 0,    0,
#'                    0.3,  0.2,  0.2,  0.3,  0,    0   ), nrow=3, byrow=TRUE)
#' 
#' test3 <- matrix(c( 0.1,  0.7,  0.2,  0,    0,    0,
#'                    0,    0,    0,    1.0,  0,    0,
#'                    0.5,  0.1,  0,    0.4,  0,    0   ), nrow=3, byrow=TRUE)
#' 
#' fromSig <- list(sig1=test1, sig2=test2, sig3=test3)
#' 
#' ### compute distances of the test signature to the converted
#' ### Alexandrov signatures from COSMIC
#' mapSignatureSets(fromSig, signAlex2Shi, method="frobenius", unique=TRUE)
#' 
#' @export mapSignatureSets
mapSignatureSets <- function(fromSignatures, toSignatures,
                             method="euclidean", unique=FALSE) {
    
    if (!isSignatureSet(fromSignatures)) {
        stop("Parameter 'fromSignatures' must be a set (list) of signatures!")
    }

    if (!isSignatureSet(toSignatures)) {
        stop("Parameter 'toSignatures' must be a set (list) of signatures!")
    }

    # from and to must be of the same format
    if (!sameSignatureFormat(fromSignatures, toSignatures)) {
        stop(paste("fromSignatures and toSignatures must be of the same",
                   "type and format!"))
    }

    # if we require a unique mapping, the number of fromSignatures must not
    # exceed the number of toSignatures
    if (unique & (length(fromSignatures) > length(toSignatures))) {
        stop(paste("for a unique mapping the number of fromSignatures",
                   "must not exceed the number of toSignatures!"))
    }

    # make sure we have signature names
    if (is.null(names(fromSignatures))) {
        names(fromSignatures) <- paste0("sign_", c(seq_along(fromSignatures)))
    }
    if (is.null(names(toSignatures))) {
        names(toSignatures) <- paste0("sign_", c(seq_along(toSignatures)))
    }

    # keep the names of the fromSignatures, but remove them from the object
    fromSigNames <- names(fromSignatures)
    names(fromSignatures) <- NULL

    
    if (!unique) {
        
        # no unique mapping is required; for each fromSignature
        # simply choose the most similar one from toSignatures
        mapping <- vapply(fromSignatures, function(from) {
            d <- determineSignatureDistances(from, toSignatures, method=method)
            names(sort(d))[1] # sort by distance, take name of first
        }, FUN.VALUE=character(1) )

        names(mapping) <- fromSigNames
    } else {
        
        # we want a unique mapping, iteratively identify the best 
        # from->to pair (shortest distance), then remove both to find the
        # next best ...
        mapping <- c()

        names(fromSignatures) <- fromSigNames
        
        while (length(fromSignatures) > 0) {

            if (length(fromSignatures) == 1 &&
                length(toSignatures) == 1) {

                # the remaining signatures are 1 for each set; this mapping
                # is clear ...
                fromS <- names(fromSignatures)
                toS <- names(toSignatures)

                # keep mapping
                names(toS) <- fromS
                mapping <- c(mapping, toS)

                # no need to continue
                break
            }

            # there are either multiple "from" or multiple "to" signatures
            # need to continue to find the minimum distance
            
            distsMatrix <- vapply(fromSignatures, function(from) {
                d <- determineSignatureDistances(from, toSignatures,
                                                 method=method)
            }, FUN.VALUE=numeric(length(toSignatures)) )

            # find minimum distance
            thisMin <- which(distsMatrix==min(distsMatrix), arr.ind=TRUE)

            fromS <- colnames(distsMatrix)[thisMin[1,"col"]]
            toS <- rownames(distsMatrix)[thisMin[1,"row"]]

            # keep mapping
            names(toS) <- fromS
            mapping <- c(mapping, toS)

            # remove the two signatures from the respective sets
            fromSignatures[thisMin[1,"col"]] <- NULL
            toSignatures[thisMin[1,"row"]] <- NULL
        }

        # reorder mapping according to fromSignature input
        mapping <- mapping[fromSigNames]
    }

    return(mapping)
}

