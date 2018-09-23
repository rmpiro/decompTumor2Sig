
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' slot (numeric value) from a \code{pmsignature::MutationFeatureData} or
#' from a \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' slot (string) from a \code{pmsignature::EstimatedParameters} or a
#' \code{pmsignature::MutationFeatureData} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
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
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSigType <- function(x) {
    stopifnot(.hasSlot(x, "type"))
    slot(x, "type")
}

#' getSigFeatDist (internal function)
#'
#' getSigFeatDist`()` extracts the content of the
#' \code{signatureFeatureDistribution} slot from a
#' \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the slot data.
#' @return Slot data.
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
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSigFeatDist <- function(x) {
    stopifnot(.hasSlot(x, "signatureFeatureDistribution"))
    slot(x, "signatureFeatureDistribution")
}


#' getSigFromEstParam (internal function)
#'
#' `getSigFromEstParam()` extracts a specific signature from the
#' \code{signatureFeatureDistribution} slot of a
#' \code{pmsignature::EstimatedParameters} object.
#' 
#' @param x Object from which to get the signature.
#' @param signum Number (1..N) of the signature to be extracted.
#' @return Signature.
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
#' @importFrom methods slot .hasSlot 
#' @keywords internal
getSigFromEstParam <- function(x, signum) {
    stopifnot(is.numeric(signum) & signum <= dim(getSigFeatDist(x))[1])
    sig <- getSigFeatDist(x)[signum,,]
}


### helper functions for checking consistency/compatibility of objects

#' isMutFeatDataObject (internal function)
#'
#' `isMutFeatDataObject()` checks whether an object is compatible to a
#' \code{pmsignature::MutationFeatureData} object (version 0.3.0), i.e.,
#' whether it contains the same info in the same format.
#' 
#' @param x Object for which to verify compatibility.
#' @return Logical value (true or false).
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
#' @importFrom methods slot .hasSlot 
#' @keywords internal
isMutFeatDataObject <- function(x) {

    # all needed data there?
    if ( !.hasSlot(x, "sampleList") ||
         !.hasSlot(x, "featureVectorList") ||
         !.hasSlot(x, "countData") ||
         !.hasSlot(x, "type") ||
         !.hasSlot(x, "flankingBasesNum") ||
         !.hasSlot(x, "transcriptionDirection")
        ) {

        return(FALSE)
    }

    # check data types
    if ( !is.character(getSampleList(x)) ||
         !is.matrix(getFeatVectList(x)) ||
         !is.matrix(getCountData(x)) ||
         !is.character(getSigType(x)) ||
         !is.integer(getNumFlankingBases(x)) ||
         !is.logical(haveTrDir(x))
        ) {
        
        return(FALSE)
    }

    numGenomes <- length(getSampleList(x))
    numFeatVecRows <- nrow(getFeatVectList(x))
    numFeatVecCols <- ncol(getFeatVectList(x))
    numCountDataRows <- nrow(getCountData(x))
    numCountDataCols <- ncol(getCountData(x))
    sigType <- getSigType(x)
    numBases <- getNumFlankingBases(x)
    trDir <- haveTrDir(x)
    
    # check consistency of the data
    if ( sigType == "independent" ) {

        # Shiraishi-type signature
        if ( (numBases + as.numeric(trDir) != numFeatVecRows) ||    # bases
             (max(getCountData(x)[1,]) != numFeatVecCols) ||        # count data
             (length(unique(getCountData(x)[2,])) != numGenomes) || # count data
             (max(getFeatVectList(x)[1,]) != 6) ||                  # feat vect.
             !(numBases == 1 || max(getFeatVectList(x)[-1,]) == 4)  # feat vect.
             
            ) {

            return(FALSE)
        }

        # check consistency in dependency of transcription direction
        if (trDir) {
            # having transcription direction
            if ( max(getFeatVectList(x)[numFeatVecRows,]) != 2 
                ) {
                return(FALSE)
            }
        }

    } else if ( sigType == "full" ) {

        # Alexandrov-type signature

        numFeatures <- 6 * 4^(numBases-1) * (as.numeric(trDir)+1)
        
        if ( (numFeatures != max(getFeatVectList(x))) ||            # bases
             (max(getCountData(x)[1,]) != numFeatVecCols) ||        # count data
             (length(unique(getCountData(x)[2,])) != numGenomes) || # count data
             (numFeatVecRows != 1) ||                               # feat vect.
             (length(getFeatVectList(x)) != numFeatVecCols)         # feat vect.
            
            ) {

            return(FALSE)
        }

    } else {
        
        # unknown singature type
        return(FALSE)
    }
    
    return(TRUE)

}


#' isEstParamObject (internal function)
#'
#' `isEstParamObject()` checks whether an object is compatible to a
#' \code{pmsignature::EstimatedParameters} object (version 0.3.0), i.e.,
#' whether it contains the same info in the same format.
#' 
#' @param x Object for which to verify compatibility.
#' @return Logical value (true or false).
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
#' @importFrom methods slot .hasSlot 
#' @keywords internal
isEstParamObject <- function(x) {

    # all needed data there?
    if ( !.hasSlot(x, "signatureNum") ||
         !.hasSlot(x, "isBackGround") ||
         !.hasSlot(x, "type") ||
         !.hasSlot(x, "signatureFeatureDistribution") ||
         !.hasSlot(x, "flankingBasesNum") ||
         !.hasSlot(x, "transcriptionDirection")
        ) {

        return(FALSE)
    }

    # check data types
    if ( !is.integer(getNumSignatures(x)) ||
         !is.logical(isBackGround(x)) ||
         !is.character(getSigType(x)) ||
         !is.array(getSigFeatDist(x)) ||
         !is.integer(getNumFlankingBases(x)) ||
         !is.logical(haveTrDir(x))
        ) {
        
        return(FALSE)
    }

    numSig <- getNumSignatures(x)
    hasBackGround <- isBackGround(x)
    sigType <- getSigType(x)
    numFeatDistDim1 <- dim(getSigFeatDist(x))[1]
    numFeatDistDim2 <- dim(getSigFeatDist(x))[2]
    numFeatDistDim3 <- dim(getSigFeatDist(x))[3]
    numBases <- getNumFlankingBases(x)
    trDir <- haveTrDir(x)
   
    # check consistency of the data
    if ( sigType == "independent" ) {
        
        # Shiraishi-type signatures
        if ( (numBases + as.numeric(trDir) != numFeatDistDim2) || # bases
             (numFeatDistDim1 != numSig) ||                       # num. sig.
             (numFeatDistDim3 != 6)                               # feat. depth
            ) {

            return(FALSE)
        }
    } else if ( sigType == "full" ) {
        
        # Alexandrov-type signatures
        if ( (6 * 4^(numBases-1) * (as.numeric(trDir)+1) != numFeatDistDim3) ||
             (numFeatDistDim1 != numSig) ||                       # num. sig.
             (numFeatDistDim2 != 1)                               # feat. depth
            ) {

            return(FALSE)
        }
    } else {
        
        # unknown singature type
        return(FALSE)
    }
    
        
    return(TRUE)

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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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

