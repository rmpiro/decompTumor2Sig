
#####################
# internal function #
#####################

#' determineTypeNumBasesAndTrDir (internal function)
#'
#' For a given signature or genome representation: determine the
#' type (Shiraishi or Alexandrov), the number of bases, and whether
#' transcription-strand information is included.
#'
#' @usage determineTypeNumBasesAndTrDir(mutData)
#' @param mutData The mutation frequency data of the signature or genome.
#' @return A list object composed of: "type"=type of signature or genome;
#' "numBases"=number of bases of the sequence patterns; and
#' "trDir"=logical value indicating whether transcription-strand information
#' has been considered in the mutation frequency data.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @keywords internal
determineTypeNumBasesAndTrDir <- function(mutData) {

    haveTrDir <- NULL
    haveNumBases <- NULL
    haveType <- NULL
    
    if(is.matrix(mutData) && ncol(mutData) == 6) {
        # Shiraishi-type

        # number of bases of sequence pattern and transcription dir.
        haveNumBases <- nrow(mutData)
        haveTrDir <- FALSE
        if ( (nrow(mutData) %% 2) == 0 ) { # even number of rows,
                                           # must have transcript direction!
            haveNumBases <- haveNumBases -1
            haveTrDir <- TRUE
        }
    
        haveType <- "Shiraishi"
        
    } else if (is.numeric(mutData)) {
        # Alexandrov-type
 
        if ( ( (log10(length(mutData)/6)/log10(4)) %% 1 ) == 0 ) {
            
            # this must be without transcription direction
            haveNumBases <- log10(length(mutData)/6)/log10(4) + 1
            haveTrDir <- FALSE
            haveType <- "Alexandrov"
            
        } else if ( ( (log10(length(mutData)/(6*2))/log10(4)) %% 1 ) == 0 ) {
            
            # this must be with transcription direction
            haveNumBases <- log10(length(mutData)/(6*2))/log10(4) + 1
            haveTrDir <- TRUE
            haveType <- "Alexandrov"
        }
        
    } else{
        warning(paste("cannot determine the number of bases and presence",
                      "of transcription direction for signature/genome!\n"))
    }

    return(list("type"=haveType, "numBases"=haveNumBases, "trDir"=haveTrDir))
}
