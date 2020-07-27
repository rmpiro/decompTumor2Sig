

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
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics 
#' 20(Suppl 4):152.\cr
#' @keywords internal
convAlx2Shi <- function(x) {

    shSig <- NULL
    
    sigprop <- determineTypeNumBasesAndTrDir(x)
    
    haveTrDir <- sigprop$trDir

    basesAdj <- (sigprop$numBases-1)/2

    strLen <- 5 + 2*basesAdj + as.numeric(haveTrDir)   # examples: "A[T>G]C"
                                                       # or "AA[T>G]AC+" 

    if (is.null(names(x))) {
        # plain vector; need to name it according to Alexandrov model
        names(x) <- buildSortedAlexandrovSignaturePatternList(sigprop$numBases,
                                                              haveTrDir)
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

