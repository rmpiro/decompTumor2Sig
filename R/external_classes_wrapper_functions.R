
######################
# internal functions #
######################


### wrapper functions for reading slot data ###

#' haveTrDir (internal function)
#'
#' `haveTrDir()` extracts the content of the \code{transcriptionDirection} slot
#' (logical value) from a \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
haveTrDir <- function(x) {
    stopifnot(.hasSlot(x, "transcriptionDirection"))
    slot(x, "transcriptionDirection")
}

#' getNumFlankingBases (internal function)
#'
#' `getNumFlankingBases()` extracts the content of the \code{flankingBasesNum}
#' slot (numeric value) from a \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getNumFlankingBases <- function(x) {
    stopifnot(.hasSlot(x, "flankingBasesNum"))
    slot(x, "flankingBasesNum")
}

#' getCountData (internal function)
#'
#' `getCountData()` extracts the content of the \code{countData}
#' slot (numeric values) from a \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getCountData <- function(x) {
    stopifnot(.hasSlot(x, "countData"))
    slot(x, "countData")
}

#' getFeatVectList (internal function)
#'
#' `getFeatVectList()` extracts the content of the \code{featureVectorList}
#' slot from a \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getFeatVectList <- function(x) {
    stopifnot(.hasSlot(x, "featureVectorList"))
    slot(x, "featureVectorList")
}

#' getSampleList (internal function)
#'
#' `getSampleList()` extracts the content of the \code{sampleList}
#' slot (strings) from a \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSampleList <- function(x) {
    stopifnot(.hasSlot(x, "sampleList"))
    slot(x, "sampleList")
}

#' getNumSignatures (internal function)
#'
#' `getNumSignatures()` extracts the content of the \code{signatureNum}
#' slot (numeric value) from a \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getNumSignatures <- function(x) {
    stopifnot(.hasSlot(x, "signatureNum"))
    slot(x, "signatureNum")
}

#' isBackGround (internal function)
#'
#' `isBackGround()` extracts the content of the \code{isBackGround}
#' slot (logical value) from a \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
isBackGround <- function(x) {
    stopifnot(.hasSlot(x, "isBackGround"))
    slot(x, "isBackGround")
}

#' getSigType (internal function)
#'
#' `getSigType()` extracts the content of the \code{type}
#' slot (string) from a \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSigType <- function(x) {
    stopifnot(.hasSlot(x, "type"))
    slot(x, "type")
}

#' getSigFeatDist (internal function)
#'
#' getSigFeatDist`()` extracts the content of the
#' \code{signatureFeatureDistribution} slot (string) from a
#' \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSigFeatDist <- function(x) {
    stopifnot(.hasSlot(x, "signatureFeatureDistribution"))
    slot(x, "signatureFeatureDistribution")
}



### wrapper functions for writing slot data ###


#' setTrDir (internal function)
#'
#' `setTrDir()` serves as wrapper functions to change content of the
#' \code{transcriptionDirection} slot (logical value) of a
#' \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object of which to change slot data.
#' @param value Value to be written to the slot.
#' @return The modified object \code{x}.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods .hasSlot 
#' @keywords internal
setTrDir <- function(x, value) {
    stopifnot(.hasSlot(x, "transcriptionDirection"))
    stopifnot(is.logical(value))
    x@transcriptionDirection <- value
    x
}

#' setNumFlankingBases (internal function)
#'
#' `setNumFlankingBases()` serves as wrapper functions to
#' change content of the \code{flankingBasesNum} slot (numeric value) of a
#' \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object of which to change slot data.
#' @param value Value to be written to the slot.
#' @return The modified object \code{x}.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods .hasSlot 
#' @keywords internal
setNumFlankingBases <- function(x, value) {
    stopifnot(.hasSlot(x, "flankingBasesNum"))
    stopifnot(is.integer(value) || (value%%2 != 1) )
    x@flankingBasesNum <- value
    x
}

#' setSigType (internal function)
#'
#' `setSigType()` serves as wrapper functions to change content of the
#' \code{type} slot (string) of a
#' \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object of which to change slot data.
#' @param value Value to be written to the slot.
#' @return The modified object \code{x}.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods .hasSlot 
#' @keywords internal
setSigType <- function(x, value) {
    stopifnot(.hasSlot(x, "type"))
    stopifnot(is.character(value))
    x@type <- value
    x
}

#' setSigFeatDist (internal function)
#'
#' `setSigFeatDist()` serves as wrapper functions to change content of the
#' \code{signatureFeatureDistribution} slot (numeric values) of a
#' \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object of which to change slot data.
#' @param value Value to be written to the slot.
#' @return The modified object \code{x}.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @importFrom methods .hasSlot 
#' @keywords internal
setSigFeatDist <- function(x, value) {
    stopifnot(.hasSlot(x, "signatureFeatureDistribution"))
    stopifnot(is.numeric(value))
    x@signatureFeatureDistribution <- value
    x
}

