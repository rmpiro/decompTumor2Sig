#' Downgrade Shiraishi-type signatures.
#'
#' `downgradeShiraishiSignatures()` downgrades/trims signatures of the
#' Shiraishi type by discarding flanking bases (reducing the length of the
#' sequence pattern) and/or the transcription direction. The downgrade doesn't
#' pose a problem because the flanking bases and the transcription direction
#' are considered as independent features according to the Shiraishi model of
#' mutational signatures.
#'
#' @usage downgradeShiraishiSignatures(signatures, numBases=NULL,
#' removeTrDir=FALSE)
#' @param signatures (Mandatory) A list of Shiraishi signatures that need to be
#' downgraded/trimmed.
#' @param numBases (Conditionally optional) The total number of bases
#' (mutated base plus flanking bases around the mutated base) that should
#' be kept. All further flanking bases farther away from the mutated bases
#' are dropped. If specified, \code{numBases} must be odd and smaller than
#' the current number of bases of the \code{signatures}. If \code{NULL}, no
#' flanking bases will be dropped. At least one of \code{numBases} or
#' \code{removeTrDir} must be specified. 
#' @param removeTrDir (Conditionally optional) Logical value that specifies
#' whether information on the transcript direction should be dropped (if
#' present at all). At least one of \code{numBases} or \code{removeTrDir}
#' must be specified.
#' @return A list of Shiraishi signatures that have been accordingly downgraded.
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
#' @seealso \code{\link{decompTumor2Sig}}
#' @examples
#' 
#' ### Load 15 Shiraishi signatures obtained from 435 tumor genomes from
#' ### Alexandrov et al. (number of bases: 5, transcription direction: yes)
#' sfile <- system.file("extdata",
#'          "Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata", 
#'          package="decompTumor2Sig")
#' load(sfile)
#' 
#' ### downgrade the signatures to include only 3 bases and drop the
#' ### transcription direction 
#' downgradeShiraishiSignatures(signatures, numBases=3, removeTrDir=TRUE)
#' 
#' @export downgradeShiraishiSignatures
downgradeShiraishiSignatures <- function(signatures, numBases=NULL,
                                         removeTrDir=FALSE) {

    if (!isShiraishiSet(signatures)) {
        stop(paste("Parameter 'signatures' must be a set (list) of Shiraishi",
                   "signatures!"))
    }

    if (is.null(removeTrDir)) {
        removeTrDir <- FALSE
    }
    if(!is.logical(removeTrDir)) {
        stop("Value of removeTrDir must be logical!")
    }

    if (is.null(numBases) & !removeTrDir) {
        stop("At least one of numBases and removeTrDir must be specified!")
    }

    if (!is.null(numBases)) {
        if (!is.numeric(numBases)
            || (numBases!=round(numBases)) || ((numBases%%2)==0) ) {
            stop("Value of numBases must be an odd integer!")
        }
    }

    newsigs <- lapply(signatures, function(sig) {
        # first, check current format of the signatures

        # if we have an even number of lines, the last one must be for the
        # transcription direction
        haveTrDir <- !(nrow(sig)%%2)

        # number of bases in the sequence pattern?
        haveBases <- nrow(sig) - as.numeric(haveTrDir)

        if(!is.null(numBases) && (numBases >= haveBases)) {
            stop(paste("Value of numBases must be smaller than the",
                       "signatures' current number of bases!"))
        }

        # determine which indices to remove
        removeRows <- c()

        if(!is.null(numBases)) {
            # upstream
            removeRows <- c(removeRows,
                            seq_len((haveBases-numBases)/2)+1)
                            #c(1:((haveBases-numBases)/2))+1)
            # downstream
            removeRows <- c(removeRows,
                            seq(from=((haveBases-(haveBases-numBases)/2)+1),
                                to=haveBases))
                           #c(((haveBases-(haveBases-numBases)/2)+1):haveBases))
        }

        if(removeTrDir & haveTrDir) {
            removeRows <- c(removeRows, nrow(sig))
        }

        sig <- sig[-removeRows,]

        if (is.vector(sig)) {
            sig <- matrix(sig, nrow=1)
        }

        setNames4ShiraishiTable(sig)
    })

    return(newsigs)
}
