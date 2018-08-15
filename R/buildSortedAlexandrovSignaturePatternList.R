
#####################
# internal function #
#####################

#' buildSortedAlexandrovSignaturePatternList (internal function)
#'
#' Build and sort the pattern list of an Alexandrov signature in the following
#' way: First according to the base change: C>A, C>G, C>T, T>A, T>C, T>G,
#' then within these categories according to the flanking bases: A, C, G, T.
#'
#' @usage buildSortedAlexandrovSignaturePatternList(numBases, trDir)
#' @param numBases Number of bases for the sequence pattern (odd integer).
#' @param trDir Logical: use transcription-strand information?
#' @return A sorted list of mutation features (e.g., triplets with base change).
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @keywords internal
buildSortedAlexandrovSignaturePatternList <- function(numBases, trDir) {
    seqpatterns <- NULL
    for (ii in seq((numBases-1),1)) {
        seqpatterns <- paste0(seqpatterns, gl(4, 4^(ii-1), 4^(numBases-1),
                                              labels=c("A","C","G","T")))
    }
    
    # add the central base
    seqpatterns <-
        paste0(substr(seqpatterns,1,(numBases-1)/2),
               sort(rep(c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"),
                        length(seqpatterns))),
               substr(seqpatterns,(numBases-1)/2+1,numBases-1))

    # double if we consider the transcription direction ...
    if (trDir) {
        seqpatterns <- paste0(c(seqpatterns,seqpatterns),
                              c(rep("+", length(seqpatterns)),
                                rep("-", length(seqpatterns))))
    }

    return(seqpatterns)
}
