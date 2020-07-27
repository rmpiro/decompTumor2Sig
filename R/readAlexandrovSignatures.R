#' Read Alexandrov-type signatures (COSMIC format).
#'
#' `readAlexandrovSignatures()` reads a set of Alexandrov-type signatures
#' (COSMIC format) from a flat file or URL. Signatures must be specified in the
#' tab-separated format used by the COSMIC website for signatures version 2
#' (March 2015) or the comma-separated format used for signatures version 3
#' (May 2019); see details below or\cr
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures_v2} ->
#' "Download signatures" for version 2, or\cr
#' \url{https://cancer.sanger.ac.uk/cosmic/signatures/SBS/} and
#' \url{https://www.synapse.org/#!Synapse:syn12009743} for version 3.
#'
#' For version 3, only Single Base Substitution (SBS) signatures can be used.
#'
#' COSMIC format for Alexandrov signatures, version 2:
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
#' COSMIC/Synapse format for Alexandrov signatures, version 3:
#'
#' \tabular{llllll}{
#' Type,SubType,SBS1,SBS2,SBS3,SBS4,SBS5,SBS6, ...\cr
#' C>A,ACA,8.86E-04,5.80E-07,2.08E-02,4.22E-02,1.20E-02,4.25E-04, ...\cr
#' C>A,ACC,2.28E-03,1.48E-04,1.65E-02,3.33E-02,9.44E-03,5.24E-04, ...\cr
#' C>A,ACG,1.77E-04,5.23E-05,1.75E-03,1.56E-02,1.85E-03,5.20E-05, ...\cr
#' C>A,ACT,1.28E-03,9.78E-05,1.22E-02,2.95E-02,6.61E-03,1.80E-04, ...\cr
#' [...]\cr
#' T>G,TTG,5.83E-04,9.54E-05,8.05E-03,2.32E-03,6.94E-03,3.24E-04, ...\cr
#' T>G,TTT,2.23E-16,2.23E-16,1.05E-02,5.68E-04,1.35E-02,1.01E-03, ...\cr
#' }
#' 
#' @usage
#' readAlexandrovSignatures(file)
#' @param file (Mandatory) Can be a file name or an URL for download.
#' Default (COSMIC):
#' "https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
#' @return A list of Alexandrov signatures that can be used for
#' \code{decomposeTumorGenomes}. 
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
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
    function(file=paste0("https://cancer.sanger.ac.uk/cancergenome/assets/",
                         "signatures_probabilities.txt")) {

    # read a set of Alexandrov signatures from a tab-separated flat file in the
    # format provided by COSMIC Mutational Signatures.
    # Version 2 signatures (March 2015) are available at: 
    #   http://cancer.sanger.ac.uk/cancergenome/assets/
    #                   signatures_probabilities.txt (this is taken by DEFAULT)
    # Version 3 signatures (May 2019) are available at:
    #   https://www.synapse.org/#!Synapse:syn12009743
    #   (only SBS signatures, i.e., for single base substitution)
        
    if (!is.character(file)) {
        stop("Parameter 'file' must be a filename or URL!")
    }

    # read all data in one table

    # first, try tab as column separator/delimiter (used for version 2):
    sigmatrix <- as.matrix(read.table(file, header=TRUE,
                                      row.names=NULL, sep="\t"))

    if (ncol(sigmatrix) < 3) {
        # cannot contain signature data; at least 2 columns of annotation
        # try comma as column separator/delimiter (used for version 3):
        sigmatrix <- as.matrix(read.table(file, header=TRUE,
                                          row.names=NULL, sep=","))
    }

    # verify format of first column (e.g., "C>A")
    if (length(grep("^[CT]>[ACGT]$", sigmatrix[,1])) != nrow(sigmatrix)) {
        stop(paste("Wrong file format. Expected SNV annotation for",
                   "pyrimidines in first column, e.g. 'C>A'."))
    }
    
    
    # verify format of second column (e.g., "ACA")
    if (length(grep("^[ACGT]{3}$", sigmatrix[,2])) != nrow(sigmatrix)) {
        stop(paste("Wrong file format. Expected mutated triplet",
                   "in second column, e.g. 'ACA'."))
    }
    
    
    # check whether we have the third annotation column (e.g., "A[C>A]A")
    if (length(grep("^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$", sigmatrix[,3])) ==
        nrow(sigmatrix)) {

        # we have the third column, use it as element names for the vector
        sigVecNames <- sigmatrix[,3]
        firstSigCol <- 4
    } else if (!anyNA(suppressWarnings(as.numeric(sigmatrix[,3])))) {
        # there are only numeric values in the third column; signature!

        # need to construct the mutation types from the first two columns
        sigVecNames <- apply(sigmatrix[,seq_len(2)], 1,
                             function(x) {
                                 y = unlist(strsplit(x[2], ""));
                                 paste0(y[1],"[",x[1],"]",y[3])
                             })
        firstSigCol <- 3
    } else {
        stop(paste("Wrong file format. Expected either a signature or",
                   "mutation type in third column, e.g. 'A[C>A]A'."))
    }

    sigList <- list()
    sigNames <- c()

    for (colId in seq(firstSigCol, ncol(sigmatrix))) {
        if (!all(is.na(sigmatrix[,colId]))) { # signature is defined

            sigVec <- as.numeric(sigmatrix[,colId])
            names(sigVec) <- sigVecNames

            # make sure this is normalized (doesn't hold for version 3
            # where the sum of probabilities is somtimes minimally different
            # from 1 ...
            sigVec <- sigVec/sum(sigVec)

            # make sure we sort the vector correctly (the version 2 file has
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


