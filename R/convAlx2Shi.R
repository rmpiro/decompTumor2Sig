

#####################
# internal function #
#####################

#' convAlx2Shi (internal function)
#'
#' Convert a single Alexandrov signature to a Shiraishi signature.
#'
#' @usage convAlx2Shi(x)
#' @param x The Alexandrov signature (mutation pattern vector) to be converted.
#' @return A Shiraishi signature (mutation pattern matrix).
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @keywords internal
convAlx2Shi <- function(x) {

    strLen <- unique(nchar(names(x)))   # example: A[T>G]C

    haveTrDir <- FALSE
    if(substr(names(x)[1], strLen, strLen) %in% c("+","-")) {
        # do we have info on transcription direction?
        haveTrDir <- TRUE
    }

    basesAdj <- (strLen-5-as.numeric(haveTrDir))/2

    # first row: base change frequencies
    shSig <- as.vector(tapply(x, substr(names(x),basesAdj+2,basesAdj+4), sum))

    # upstream base frequencies
    for (ii in seq(0,(basesAdj-1))) {
        shSig <- rbind(shSig,
            c(as.vector(tapply(x, substr(names(x),1+ii,1+ii), sum)), rep(0,2))
                       )
    }

    # downstream base frequencies
    for (ii in seq(0,(basesAdj-1))) {
        shSig <- rbind(shSig,
                       c(as.vector(tapply(x, substr(names(x),
                                                    strLen-ii
                                                        -as.numeric(haveTrDir),
                                                    strLen-ii
                                                        -as.numeric(haveTrDir)),
                                          sum)), rep(0,2))
                       )
    }

    if(haveTrDir) {
        shSig <- rbind(shSig,
                       c(as.vector(tapply(x, substr(names(x),strLen,strLen),
                                          sum)[c("+","-")]), rep(0,4))
                       )
    }

    rownames(shSig) <- NULL

    return(shSig)
}

