#' Read Alexandrov-type signatures (COSMIC format).
#'
#' `readAlexandrovSignatures()` reads a set of Alexandrov-type signatures
#' (COSMIC format) from a flat file or URL. Signatures must be specified in the
#' tab-separated format used by the COSMIC website; see details below or\cr
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures} -> "Download signatures". 
#'
#' COSMIC format for Alexandrov signatures:
#'
#' \tabular{llllll}{
#' Subst. \tab Trinucleotide \tab Mutation Type \tab Signature 1 \tab
#' Signature 2 \tab...\cr
#' C>A \tab ACA \tab A[C>A]A \tab 0.0110983262 \tab 0.0006827082 \tab...\cr
#' C>A \tab ACC \tab A[C>A]C \tab 0.0091493407 \tab 0.0006191072 \tab...\cr
#' C>A \tab ACG \tab A[C>A]G \tab 0.0014900705 \tab 0.0000992790 \tab...\cr
#' C>A \tab ACT \tab A[C>A]T \tab 0.0062338852 \tab 0.0003238914 \tab...\cr
#' [...]\cr
#' T>G \tab TTG \tab T[T>G]G \tab 0.0020310769 \tab 0.0002066152 \tab...\cr
#' T>G \tab TTT \tab T[T>G]T \tab 0.0040301281 \tab 0.0000235982 \tab...\cr
#' }
#' 
#' @usage
#' readAlexandrovSignatures(file)
#' @param file (Mandatory) Can be a file name or an URL for download.
#' Default (COSMIC):
#' "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
#' @return A list of Alexandrov signatures that can be used for
#' \code{decomposeTumorGenomes}. 
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
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{readShiraishiSignatures}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signatures <- readAlexandrovSignatures()
#' 
#' @importFrom utils read.table
#' @export readAlexandrovSignatures
readAlexandrovSignatures <-
    function(file=paste0("http://cancer.sanger.ac.uk/cancergenome/assets/",
                 "signatures_probabilities.txt")) {

    # read a set of Alexandrov signatures from a tab-separated flat file in the
    # format provided by COSMIC Mutational Signatures, available at: 
    # http://cancer.sanger.ac.uk/cancergenome/assets/
    #                   signatures_probabilities.txt (this is taken by default)
 
    if (!is.character(file)) {
        stop("Parameter 'file' must be a filename or URL!")
    }

    # read all data in one table
    sigmatrix <- as.matrix(read.table(file, header=TRUE,
                                      row.names=NULL, sep="\t"))

    sigList <- list()
    sigNames <- c()

    for (colId in seq(4,ncol(sigmatrix))) {
        if (!all(is.na(sigmatrix[,colId]))) { # signature is defined
            sigVec <- as.numeric(sigmatrix[,colId])
            names(sigVec) <- sigmatrix[,3]

            # make sure we sort the vector correctly (the official file has
            # another sorting): first check number of bases and presence of
            # transcription direction
            sigFeatures <- determineTypeNumBasesAndTrDir(sigVec)
            if (is.null(sigFeatures$type) || sigFeatures$type != "Alexandrov") {
                stop("Wrong number of patterns for an Alexandrov signature!")
            }
            # reorder the signature vector so that base changes stay together
            sigVec <-
                sigVec[buildSortedAlexandrovSignaturePatternList(
                    numBases=sigFeatures$numBases,
                    trDir=sigFeatures$trDir)]
            
            # the vector is fine now, save in list of multiple signatures
            sigList[[length(sigList)+1]] <- sigVec
            sigNames <- c(sigNames, colnames(sigmatrix)[colId])
        }
    }

    names(sigList) <- sigNames

    return(sigList)
}


