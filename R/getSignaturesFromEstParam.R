#' Get signatures from an \code{EstimatedParameters} object.
#'
#' `getSignaturesFromEstParam()` takes an \code{EstimatedParameters} object
#' (signatures data) as computed by the '\code{pmsignature}' package (by
#' \code{pmsignature::getPMSignature}; version 0.3.0) and extracts the
#' signature information. The signatures can then be passed to
#' \code{decomposeTumorGenomes}. 
#'
#' @usage getSignaturesFromEstParam(Param)
#' @param Param (Mandatory) A \code{pmsignature::EstimatedParameters} object
#' as those produced by the de novo signature construction method
#' \code{pmsignature::getPMSignature}.
#' @return A list of Shiraishi signatures, one object per signature. Please
#' see \code{readShiraishiSignatures} or the \code{decompTumor2Sig} vignette
#' for more information on the format of Shiraishi signatures.
#' @author Rosario M. Piro and Sandra Krueger\cr Freie Universitaet Berlin\cr
#' Maintainer: Rosario M. Piro\cr E-Mail: <rmpiro@@gmail.com> or
#' <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{readShiraishiSignatures}}
#' @examples
#' 
#' ### load signatures for breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) in the format produced by
#' ### pmsignature (PMID: 26630308)
#' pmsigdata <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-pmsignature-Param.Rdata", 
#'          package="decompTumor2Sig")
#' load(pmsigdata)
#' 
#' ### extract the signatures from the pmsignature Param object
#' signatures <- getSignaturesFromEstParam(Param)
#' 
#' @importFrom methods is
#' @export getSignaturesFromEstParam
getSignaturesFromEstParam <- function(Param) {
    # convert the signatures from an "EstimatedParameters" (Param),
    # as obtained from pmsignature, to a list of signature matrices or vectors

    if (!isEstParamObject(Param)) {
        stop(paste("Param must be compatible to an object of type",
                   "EstimatedParameters (as produced by pmsignature's ",
                   "function getPMSignature; version 0.3.0)"))
    }

    sigList <- list()

    # number of signatures, etc.
    numSigs <- getNumSignatures(Param)
    numBases <- getNumFlankingBases(Param)
    trDir <- haveTrDir(Param)
    
    # if one of them was background, the true number of signatures is one less!
    numSigs <- numSigs - as.numeric(isBackGround(Param))

    for (sig in seq_len(numSigs)) {
        sigList[[length(sigList)+1]] <-
            getSigFromEstParam(Param, sig)

        # set row and/or column names
        if (getSigType(Param) == "independent") {
            # Shiraishi
            if (is.vector(sigList[[length(sigList)]])) {
                # this can happen when numBases=1 and trDir=FALSE; would be vector, we need table!
                sigList[[length(sigList)]] <- matrix(sigList[[length(sigList)]], nrow=1)
            }
            
            sigList[[length(sigList)]] <-
                setNames4ShiraishiTable(sigList[[length(sigList)]])
        } else {
            # Alexandrov
            names(sigList[[length(sigList)]]) <-
                buildSortedAlexandrovSignaturePatternList(numBases=numBases,
                                                          trDir=trDir)
        }
    }

    return(sigList)
}


