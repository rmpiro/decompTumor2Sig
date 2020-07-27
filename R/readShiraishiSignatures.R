#' Read a set of Shiraishi signatures.
#'
#' `readShiraishiSignatures()` reads one or more Shiraishi-type signatures
#' from flat files (one file per signature). The signatures must be specified
#' as matrices without headers and row names (see details below).
#'
#' Format (see Shiraishi et al. PLoS Genetics 11(12):e1005657, 2015):
#' 
#' First line: Frequencies of the base changes C>A, C>G, C>T, T>A, T>C,
#' and T>G
#' 
#' Following 2k lines (for k up- and downstream flanking bases):
#' Frequencies of the bases A, C, G, and T, followed by two 0 values
#' 
#' Final line (only if transcription direction is considered):
#' Frequencies of occurrences on the transcription strand, and on the
#' opposite strand, followed by four 0 values.
#' 
#' Example:
#' \tabular{llllll}{
#' 1.8874e-14 \tab 0.10974 \tab 0.045918 \tab
#' 0.11308 \tab 0.07429 \tab 0.65697\cr
#' 3.8079e-01 \tab 0.12215 \tab 0.191456 \tab
#' 0.30561 \tab 0.00000 \tab 0.00000\cr
#' 1.5311e-01 \tab 0.34214 \tab 0.179774 \tab
#' 0.32497 \tab 0.00000 \tab 0.00000\cr
#' 1.2378e-01 \tab 0.10243 \tab 0.163461 \tab
#' 0.61032 \tab 0.00000 \tab 0.00000\cr
#' 3.4891e-01 \tab 0.15346 \tab 0.156687 \tab
#' 0.34094 \tab 0.00000 \tab 0.00000\cr
#' 5.6435e-01 \tab 0.43565 \tab 0.000000 \tab
#' 0.00000 \tab 0.00000 \tab 0.00000 \cr
#' }
#'
#' @usage readShiraishiSignatures(files)
#' @param files (Mandatory) Can be a single file name, a vector of file names,
#' or a list of file names.
#' @return A list of Shiraishi signatures that can be used for
#' \code{decomposeTumorGenomes}.
#' @author Rosario M. Piro, Politecnico di Milano\cr
#' Sandra Krueger, Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{readAlexandrovSignatures}}\cr
#' \code{\link{getSignaturesFromEstParam}}
#' @examples
#' 
#' ### read four Shiraishi signatures for breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) from flat files
#' sigfiles <- system.file("extdata",
#'          paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",1:4,".tsv"),
#'          package="decompTumor2Sig")
#' 
#' signatures <- readShiraishiSignatures(sigfiles)
#' 
#' @importFrom utils read.table
#' @export readShiraishiSignatures
readShiraishiSignatures <- function(files) {
    # read a set of Shiraishi signatures from flat files
    # files can be a single file name of vector of file names as returned
    #   by list.files(path="", pattern="" ....)

    if (is.list(files)) {
        files <- unlist(files)
    }

    if (!is.character(files)) {
        stop(paste("Parameter 'files' must be a filename or a list/vector",
                   "of file names!"))
    }
    
    fnames <- as.vector(files)
    
    sigList <- list()

    for (fn in fnames) {
        sigmatrix <- as.matrix(read.table(fn, header=FALSE, row.names=NULL))

        sigmatrix <- setNames4ShiraishiTable(sigmatrix)

        sigList[[length(sigList)+1]] <- sigmatrix
    }

    if (length(unique(basename(fnames))) == length(fnames)) {
        # base file names (without paths) are unique, take as signature names
        names(sigList) <- basename(fnames)
    } else {
        # take the unique path/file names as signature names
        names(sigList) <- fnames
    }
    
    return(sigList)
}

