
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
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics 
#' 20(Suppl 4):152.\cr
#' @keywords internal
determineTypeNumBasesAndTrDir <- function(mutData) {

    haveTrDir <- NULL
    haveNumBases <- NULL
    haveType <- NULL
    
    if(isShiraishiSet(list(mutData))) {
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
        
    } else if (isAlexandrovSet(list(mutData))) {
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
