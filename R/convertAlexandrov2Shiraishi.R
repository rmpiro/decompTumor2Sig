#' Convert Alexandrov-type signatures to Shiraishi signatures
#'
#' `convertAlexandrov2Shiraishi()` converts a set Alexandrov signatures to the
#' Shiraishi model, summing the respective frequencies of base
#' changes, and upstream and downstream flanking bases. In most cases, the
#' resulting Shiraishi signatures don't provide information on the
#' transcription strand, as this is not part of the standard Alexandrov
#' signatures. While the conversion is mainly thought for signatures, it
#' actually works also for mutation frequency data from genomes which have
#' the same format. [Attention: this conversion entails a loss of specificity
#' and the applicability of Shiraishi signatures derived from Alexandrov
#' signatures has not been extensively explored!]
#'
#' @usage convertAlexandrov2Shiraishi(signatures)
#' @param signatures (Mandatory) A list of Alexandrov signatures with named
#' elements as produced by \code{readAlexandrovSignatures}.
#' @return A list of Shiraishi signatures that can be used for
#' \code{decomposeTumorGenomes}.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr
#' Maintainer: Rosario M. Piro\cr E-Mail: <r.piro@@fu-berlin.de> or
#' <rmpiro@@gmail.com>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{readAlexandrovSignatures}}\cr
#' \code{\link{readShiraishiSignatures}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' ### convert them to the Shiraishi model
#' signShiraishi <- convertAlexandrov2Shiraishi(signAlexandrov)
#' 
#' @export convertAlexandrov2Shiraishi
convertAlexandrov2Shiraishi <- function(signatures) {

    if (is.probability.vector(signatures)) {
        # if the user specified a single signature, make it a list
        signatures <- list(signatures)
    }

    if (!isAlexandrovSet(signatures)) {
        stop(paste("Parameter signatures must be a set (list) of Alexandrov",
                   "signatures."))
    }

    shSignatures <- lapply(signatures, convAlx2Shi)

    return(shSignatures)
}

