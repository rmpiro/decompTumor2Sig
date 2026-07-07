
#####################
# internal function #
#####################

#' convertMutTypesFromAlexandrov_v3_5 (internal function)
#'
#' Convert the mutation type naming from the "NNNN" to "N[N>N]N" (e.g.,
#' "ACGT" corresponds to "A[C>T]G"!)
#' 
#' @usage convertMutTypesFromAlexandrov_v3_5(mutNames)
#' @param mutNames A vector of mutation type names of format "NNNN"
#' @return A vector of mutation type names of format "N[N>N]N"
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
convertMutTypesFromAlexandrov_v3_5 <- function(mutNames) {

    if (!is.character(mutNames) || length(grep("^[ACGT][ACGT][ACGT][ACGT]$",
                                               mutNames)) != length(mutNames)
        ) {
        stop (paste("Couldn't convert mutation-type names from \"NNNN\" to",
                    "\"N[N>N]N\"."))
    }

    leftBase <- substr(mutNames, 1, 1)
    mutBase <- substr(mutNames, 2, 2)
    rightBase <- substr(mutNames, 3, 3)
    altBase <- substr(mutNames, 4, 4)
    
    mutNames <- paste0(leftBase,"[",mutBase,">",altBase,"]",rightBase)

    return(mutNames)
}


