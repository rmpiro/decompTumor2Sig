
#####################
# internal function #
#####################

#' QPforSig (internal function)
#'
#' Perform quadratic programming for signatures to determine exposures.
#'
#' @usage QPforSig(counts, signatures, constrainToMaxContribution=FALSE,
#' tolerance=0.1)
#' @param counts The genome's mutation frequencies (either in Alexandrov or
#' Shiraishi format).
#' @param signatures The signatures to be used for decomposition. Must be in
#' the same format as the genome's \code{counts}.
#' @param constrainToMaxContribution (Optional) [Note: this is experimental
#' and is usually not needed!] If \code{TRUE}, the maximum contribution that
#' can be attributed to a signature will be constraint by the variant feature
#' counts (e.g., specific flanking bases) observed in the individual tumor
#' genome. If, for example, 30\% of all observed variants have a specific
#' feature and 60\% of the variants produced by a mutational process/signature
#' will manifest the feature, then the signature can have contributed up to
#' 0.3/0.6 (=0.5 or 50\%) of the observed variants. The lowest possible
#' contribution over all signature features will be taken as the allowed
#' maximum contribution of the signature. This allowed maximum will
#' additionally be increased by the value specified as \code{tolerance}
#' (see below). For the illustrated example and \code{tolerance}=0.1 a
#' contribution of up to 0.5+0.1 = 0.6 (or 60\%) of the signature would be
#' allowed. 
#' @param tolerance (Optional) If \code{constrainToMaxContribution} is
#' \code{TRUE}, the maximum contribution computed for a signature is increased
#' by this value (see above). If the parameter \code{constrainToMaxContribution}
#' is \code{FALSE}, the tolerance value is ignored. Default: 0.1.
#' @return The decomposition in form of an exposure vector (same order as
#' \code{signatures}).
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom Matrix nearPD
#' @keywords internal
QPforSig <- function(counts, signatures, constrainToMaxContribution=FALSE,
                     tolerance=0.1) {

    if (!isSignatureSet(signatures)) {
        stop("Signatures must be a set (list) for function QPforSig")
    }

    if (constrainToMaxContribution && (tolerance < 0 ||  tolerance > 1)) {
        stop(paste("Tolerance must be between 0 and 1 when constraining the",
                   "maximum contribution of signatures in function QPforSig"))
    }

    # if necessary, convert signatures from matrices to vectors
    if (isShiraishiSet(signatures)) {
        
        signatures <-
            lapply(signatures, function(x){as.vector(t(as.matrix(x)))})
    }

    # if necessary, convert counts from matrix to vector
    if(isShiraishiSet(list(counts))) {
        counts <- as.vector(t(as.matrix(counts)))
    }

    # the signatures should be stored in a matrix, one row is one signature
    sigMa <- do.call(rbind,signatures) 

    ### needed to make sure that the matrix (sigMa %*% t(sigMa)) is
    ### positive definite!
    sigMaSquared.nearPD =
        nearPD((sigMa %*% t(sigMa)), corr=FALSE, keepDiag=TRUE,
               ensureSymmetry=TRUE)$mat

    #Rinverse <- backsolve(chol(sigMa %*% t(sigMa)),diag(dim(sigMa)[1]))
    Rinverse <- backsolve(chol(sigMaSquared.nearPD),diag(dim(sigMa)[1]))

    if (constrainToMaxContribution) {
        # We want to constrain the maximum contribution each signature can
        # have (given some tolerance):
        # For each row in sigMa (each signature), we divide the counts from
        # the genome by the corresponding fractions in the signature and take
        # the minimum over the signature, the signature can't have a higher
        # contribution than this ratio, but we add a tolerance value (given
        # that the data is noisy!)
        # Example: if in the genome 30% of variants have a specific feature
        # (e.g. a specific flanking base), i.e., the counts score is 0.3, and
        # 60% (0.6) of the variants produced by a process/signature have this 
        # feature, then the signature can have contributed up to 0.3/0.6 = 0.5
        # (50%) of the genome's variants. 
        # With the default tolerance of 0.1, we would allow the quadratic
        # programming approach to assign up to 0.5+0.1 (60%) of the variants
        # to this signature

        maxContributions <-
            apply(sigMa, 1, function(x) {
                min((counts/x)[!is.na(counts/x)]) + tolerance })

        # limit maximum contribution to 100%
        maxContributions[maxContributions>1] <- 1


        # for the additional constraints, we use negative values because we 
        # require contribution/weight w <= w_max, which we can implement 
        # as -w >= -w_max (because quadprog uses Ax >= b, so we have to
        # use >= instead of <=!)
        Amat <- cbind(rep(1,dim(Rinverse)[1]),
                      diag(dim(Rinverse)[1]), -diag(length(maxContributions)))
        bvec <- c(1, rep(0,dim(Rinverse)[1]), -maxContributions)
    }
    else { # !constrainToMaxContribution
        Amat <- cbind(rep(1,dim(Rinverse)[1]), diag(dim(Rinverse)[1]))
        bvec <- c(1, rep(0,dim(Rinverse)[1]))
    }

    dvec <- (counts) %*% t(sigMa)
    resQP <- quadprog::solve.QP(Dmat=Rinverse, dvec=dvec, Amat=Amat,
                                bvec=bvec, meq=1, factorized=TRUE)
    resQP$solution[resQP$solution<0]=0  # there can be some values below 0
                                        # which should actually be 0 (very 
                                        # close to 0 in any case, e.g., -10^-15)

    return(resQP$solution)
}

