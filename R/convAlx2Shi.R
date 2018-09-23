

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
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @keywords internal
convAlx2Shi <- function(x) {

    shSig <- NULL
    
    sigprop <- determineTypeNumBasesAndTrDir(x)
    
    haveTrDir <- sigprop$trDir

    basesAdj <- (sigprop$numBases-1)/2

    strLen <- 5 + 2*basesAdj + as.numeric(haveTrDir)   # examples: "A[T>G]C" or "AA[T>G]AC+" 

    if (is.null(names(x))) {
        # plain vector; need to name it according to Alexandrov model
        names(x) <- buildSortedAlexandrovSignaturePatternList(sigprop$numBases, haveTrDir)
    }
    
    # first row: base change frequencies
    shSig <- rbind(shSig, as.vector(tapply(x, substr(names(x),
                                                     basesAdj+2,
                                                     basesAdj+4),
                                           sum)))

    if (basesAdj > 0) {
        # upstream base frequencies
        for (ii in seq(0,(basesAdj-1))) {
            shSig <- rbind(shSig,
                           c(as.vector(tapply(x, substr(names(x),1+ii,1+ii),
                                              sum)), rep(0,2))
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
    }

    if(haveTrDir) {
        shSig <- rbind(shSig,
                       c(as.vector(tapply(x, substr(names(x),strLen,strLen),
                                          sum)[c("+","-")]), rep(0,4))
                       )
    }

    setNames4ShiraishiTable(shSig)
}

