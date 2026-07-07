#' Read Alexandrov-type signatures (COSMIC format).
#'
#' `readAlexandrovSignatures()` reads a set of Alexandrov-type signatures
#' (COSMIC format) from a flat file or URL. Signatures must be specified in the
#' tab-separated format used by the COSMIC website for signatures version 2
#' (March 2015), the comma-separated format used for signatures version 3
#' (May 2019), the Microsoft Excel 2007+ sheet used for version 3.1, the
#' tab-sperated format used for version 3.2 onwards, or the format used
#' by the CRAN package `cosmicsig` (see Details below).
#' Excel sheets cannot be read from an URL and must be downloaded first.
#'
#' Alternatively, instead of a file name `readAlexandrovSignatures()` can also
#' accept a numeric matrix object with signatures in columns (signature
#' names taken from column names) and mutation types in rows (row names
#' must correspond to the mutation-type naming of signature version 3.2
#' or of package cosmicsig, see below). This is compatible with the S3-class
#' `SBS96Catalog` format as provided by the CRAN package `cosmicsig`. 
#' 
#' For details on the accepted signature formats, see below or\cr
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures_v2} ->
#' "Download signatures" for version 2,\cr
#' \url{https://www.synapse.org/#!Synapse:syn12009743} for version 3,\cr
#' \url{https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx}
#' for version 3.1, \cr
#' and \url{https://cancer.sanger.ac.uk/signatures/} for version 3.2 onwards. 
#' For versions 3 and above, only Single Base Substitution (SBS)
#' signatures can be used.
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
#' COSMIC/Synapse format for Alexandrov signatures, version 3 and 3.1:
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
#' Version 3.1 has assentially the same format as version 3, but is distributed
#' as an Excel spread sheet.
#' 
#' COSMIC/Synapse format for Alexandrov signatures, version 3.2 onwards:
#'
#' \tabular{llll}{
#' Type \tab SBS1 \tab SBS2 \tab...\cr
#' A[C>A]A \tab 0.0110983262 \tab 0.0006827082 \tab...\cr
#' A[C>A]C \tab 0.0091493407 \tab 0.0006191072 \tab...\cr
#' A[C>A]G \tab 0.0014900705 \tab 0.0000992790 \tab...\cr
#' A[C>A]T \tab 0.0062338852 \tab 0.0003238914 \tab...\cr
#' [...]\cr
#' T[T>G]G \tab 0.0020310769 \tab 0.0002066152 \tab...\cr
#' T[T>G]T \tab 0.0040301281 \tab 0.0000235982 \tab...\cr
#' }
#'
#' COSMIC format for Alexandrov signatures as provided by `cosmicsig`:
#'
#' \tabular{llll}{
#' Type \tab SBS1 \tab SBS2 \tab...\cr
#' ACAA \tab 0.0110983262 \tab 0.0006827082 \tab...\cr
#' ACCA \tab 0.0091493407 \tab 0.0006191072 \tab...\cr
#' ACGA \tab 0.0014900705 \tab 0.0000992790 \tab...\cr
#' ACTA \tab 0.0062338852 \tab 0.0003238914 \tab...\cr
#' [...]\cr
#' TTGG \tab 0.0020310769 \tab 0.0002066152 \tab...\cr
#' TTTG \tab 0.0040301281 \tab 0.0000235982 \tab...\cr
#' }
#' (Note: the first three nucleotides are the mutated triplet, the forth
#' nucleotide is the alternate base, e.g. "ACGT" corresponds to "A[C>T]G"!)
#'
#' @usage
#' readAlexandrovSignatures(file, sigmatrix)
#' @param file (Mandatory if \code{sigmatrix} not specified) Can be a file
#' name or an URL for download.
#' @param sigmatrix (Mandatory if \code{file} not specified) A numeric matrix
#' (with named rows and columns) containing mutation types in rows and
#' individual signatures in columns. Row names must use the version 3.2 or 3.5
#' naming of mutation types (see above). \code{SBS96Catalog} matrices as
#' provided by the CRAN package \code{cosmicsig} can be used. 
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
#' ### get Alexandrov signatures v 3.5 from cosmicsig
#' library(cosmicsig)
#' sigm <- COSMIC_v3.5$signature$GRCh37$SBS9
#'
#' ### exclude possible artifacts
#' sigm <- sigm[, !(colnames(sigm) %in% possible_artifacts())]
#'
#' ### read/import signatures from the cosmicsig format
#' signatures <- readAlexandrovSignatures(sigmatrix=sigm)
#' 
#' @importFrom utils read.table
#' @importFrom readxl read_excel excel_format
#' @export readAlexandrovSignatures
readAlexandrovSignatures <- function(file=NULL, sigmatrix=NULL) {

    sigVecNames <- NULL
    firstSigCol <- NULL

    # read a set of Alexandrov signatures from a tab-separated flat file in one
    # of the formats that are (or were) provided by COSMIC Mutational Signatures
    # or from a matrix (see above)

    if (is.null(file) && is.null(sigmatrix)) {
        stop("Either 'file' or 'sigmatrix' must be specified.")
    } else if (!is.null(file) && !is.null(sigmatrix)) {
        stop("Use either a file name or a signature matrix, not both.")
    }


    # If a file has been specified, try to read the signature matrix from there

    if (!is.null(file)) { # data needs to be read from file or URL
        
        if (!is.character(file)) {
            stop("Parameter 'file' must be a file name or an URL")
        }
 
        # First, try to read this as an Excel sheet (version 3.1):
        if (!is.na(readxl::excel_format(file))) {
            sigmatrix <- readxl::read_excel(file, trim_ws=TRUE)

            # convert to data frame
            sigmatrix <- as.data.frame(sigmatrix)

        } else {
            # This is NOT an Excel file, so read it as TSV or CSV table!

            sigmatrix <- tryCatch(
                # try tab as column separator/delimiter (version 2 and 3.2):
                as.matrix(read.table(file, header=TRUE,
                                     row.names=NULL, sep="\t")),
                warning = function(cond) {
                    return(NULL)
                },
                error = function(cond) {
                    return(NULL)
                }
            )

            if (is.null(sigmatrix) || ncol(sigmatrix) < 2) {
                # cannot contain signature data (at least 1 column of anno)
                # try comma as column separator/delimiter (used for version 3):
                sigmatrix <- tryCatch(
                    as.matrix(read.table(file, header=TRUE,
                                         row.names=NULL, sep=",")),
                    warning = function(cond) {
                        return(NULL)
                    },
                    error = function(cond) {
                        return(NULL)
                    }
                )
            }
        }
        
        if (is.null(sigmatrix) || ncol(sigmatrix) < 2) {
            stop(paste("Couldn't read", file, "as tabular data file."))
        }

        
        # Check that we have one of the expected formats!
        # This was read from a file, mutation types are somewhere in the first 
        # column(s)

        if (length(grep("^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$", sigmatrix[,1])) ==
            nrow(sigmatrix) || length(grep("^[ACGT][ACGT][ACGT][ACGT]$",
                                           sigmatrix[,1])) == nrow(sigmatrix)
            ) {
            # This is version 3.2 or 3.5+, there is only one annotation column
        
            sigVecNames <- sigmatrix[,1]  # take first column as mutation types
            firstSigCol <- 2              # data starts from the second column

        } else if (ncol(sigmatrix) > 2
                   && length(grep("^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$",
                                  sigmatrix[,3])) == nrow(sigmatrix)
                   ) {
            # This is verion 2, the mutation type is in the third column
        
            # we have the third column, use it as element names for the vector
            sigVecNames <- sigmatrix[,3]
            firstSigCol <- 4

        } else {
            # This might be version 3 or 3.1; get mutation type from cols 1&2

            # verify format of first column (e.g., "C>A")
            # verify format of second column (e.g., "ACA")
            if (length(grep("^[CT]>[ACGT]$", sigmatrix[,1])) != nrow(sigmatrix)
                || length(grep("^[ACGT]{3}$", sigmatrix[,2]))
                != nrow(sigmatrix)) {
                stop(paste("Wrong file format. Need mutation type in first",
                           "or third column (e.g., 'A[C>A]A'); and/or SNV",
                           "annotation for pyrimidines in first column (e.g.",
                           "'C>A') and mutated triplet in second column (e.g.",
                           "'ACA')."))
            } else {
                # need to construct the mutation types from the first two cols
                sigVecNames <- apply(sigmatrix[,seq_len(2)], 1,
                                     function(x) {
                                         y <- unlist(strsplit(x[2], ""));
                                         paste0(y[1],"[",x[1],"]",y[3])
                                     })
                firstSigCol <- 3
            }
        }

        # end reading of signature file

    } else { # not a file but a matrix has been specified, take mutation names
             # from row names

        if (!is.matrix(sigmatrix) || !is.numeric(sigmatrix)) {
            stop("Parameter 'sigmatrix' must be a numeric matrix!")
        }

        sigVecNames <- rownames(sigmatrix)
        firstSigCol <- 1
    }


    # Do the mutation-type names need conversion?
    if (length(grep("^[ACGT][ACGT][ACGT][ACGT]$",
                    sigVecNames)) == nrow(sigmatrix)) {

        # convert the mutation type naming from the "NNNN" to "N[N>N]N"
        sigVecNames <- convertMutTypesFromAlexandrov_v3_5(sigVecNames)
    }


    # Final checks: 1) do we have correct mutation-type names? need conversion?
    if (length(grep("^[ACGT]\\[[CT]>[ACGT]\\][ACGT]$", sigVecNames)) !=
        nrow(sigmatrix) ) {

        stop("Could not get readable mutation-type names for signatures!")
    }

    # Final checks: 2) do we have numeric data in the next column?
    if (anyNA(suppressWarnings(as.numeric(sigmatrix[,firstSigCol])))) {
        stop(paste0("Wrong file or matrix format. Expected a numeric ",
                    "signature in column ",firstSigCol,"."))
    }

    # Final checks: 3) do we have column names to get signatures names? 
    if (is.null(colnames(sigmatrix))) {
        stop("Could not get names for signatures (no column names present)!")
    }

    
    # We have, what we need and can extract the signatures as a list!
    
    # now extract the signatures from the table
    sigList <- list()
    sigNames <- c()

    for (colId in seq(firstSigCol, ncol(sigmatrix))) {
        if (!all(is.na(sigmatrix[,colId]))) { # signature is defined

            sigVec <- as.numeric(sigmatrix[,colId])
            names(sigVec) <- sigVecNames

            # make sure this is normalized (doesn't hold for versions 3 and 3.1
            # where the sum of probabilities is somtimes minimally different
            # from 1 ...
            while(sum(sigVec) != 1) {  # once wasn't always sufficient!
                sigVec <- sigVec/sum(sigVec)
            }

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


